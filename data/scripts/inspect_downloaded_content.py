#!/usr/bin/env python3
"""Describe complete downloaded files without joining experiments or fitting models."""
import csv
import gzip
import json
from pathlib import Path
import zipfile

import h5py
import pyarrow.parquet as pq

DATA = Path(__file__).resolve().parents[1]
status = json.loads((DATA / "download_status.json").read_text())
records = []
previous_file = DATA / "下载文件内容清单.csv"
previous = {}
if previous_file.exists():
    with previous_file.open(encoding="utf-8-sig", newline="") as f:
        previous = {row["path"]: row for row in csv.DictReader(f)}
for item in status["files"]:
    if item["status"] != "已完成":
        continue
    path = DATA / item["path"]
    cached = previous.get(item["path"])
    if cached and int(cached["bytes"]) == path.stat().st_size:
        records.append(cached)
        continue
    out = dict(dataset=item["dataset"], file=item["name"], role=item["role"], purpose=item["purpose"],
               path=item["path"], bytes=path.stat().st_size, scope=item["scope"])
    name = path.name
    if name.endswith(".csv"):
        with path.open(encoding="utf-8-sig", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = sum(1 for _ in reader)
        out.update(rows=rows, columns=len(header), fields_sample="; ".join(header[:8]), content_check="CSV完整遍历可读")
    elif name.endswith(".h5ad"):
        with h5py.File(path, "r") as f:
            shape = f["X"].attrs.get("shape")
            if shape is None and isinstance(f["X"], h5py.Dataset):
                shape = f["X"].shape
            out.update(rows=int(shape[0]), columns=int(shape[1]), fields_sample="; ".join(list(f["obs"])[:10]),
                       content_check="HDF5表达与观测字段可读；不做分析质控")
    elif name.endswith(".parquet"):
        meta = pq.read_metadata(path)
        out.update(rows=meta.num_rows, columns=meta.num_columns, fields_sample="; ".join(meta.schema.names[:10]), content_check="Parquet元数据可读")
    elif name.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            members = [x for x in z.infolist() if not x.is_dir() and not x.filename.startswith("__MACOSX")]
            out.update(archive_files=len(members), unpacked_bytes=sum(x.file_size for x in members),
                       fields_sample="; ".join(x.filename for x in members[:8]), content_check="下载时已通过ZIP CRC检查")
    elif name.endswith(".gctx.gz"):
        with gzip.open(path, "rb") as f:
            magic = f.read(8)
        out.update(fields_sample=magic.hex(), content_check="官方SHA512及压缩流头已检查；保留完整压缩矩阵")
    elif name.endswith("series_matrix.txt.gz"):
        with gzip.open(path, "rt", errors="strict") as f:
            first = f.readline().strip()
        out.update(fields_sample=first[:180], content_check="下载时已遍历gzip检查CRC；GEO文本可读")
    else:
        out["content_check"] = "来源文件已完成大小或校验码与格式检查；未进行实验数据对接"
    records.append(out)
keys = list(dict.fromkeys(k for row in records for k in row))
with (DATA / "下载文件内容清单.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(records)
print("Content descriptions saved for", len(records), "complete files; no datasets joined")
