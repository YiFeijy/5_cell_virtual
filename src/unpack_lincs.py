#!/usr/bin/env python3
"""Expand fully verified LINCS matrices; never read incomplete .part files."""
import gzip
import json
from pathlib import Path
import shutil

import h5py

DATA = Path(__file__).resolve().parents[1] / 'data'
status = json.loads((DATA / "download_status.json").read_text())
records = []
for item in status["files"]:
    if item["status"] != "已完成" or not item["name"].endswith(".gctx.gz"):
        continue
    if not item.get("official_checksum_verified"):
        raise ValueError("Official checksum has not been verified: " + item["id"])
    source = DATA / item["path"]
    dest = DATA / "unpacked" / item["dataset"] / source.name.removesuffix(".gz")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        part = dest.with_name(dest.name + ".unpacking")
        with gzip.open(source, "rb") as reader, part.open("wb") as writer:
            shutil.copyfileobj(reader, writer, length=8 * 1024 * 1024)
        part.replace(dest)
    entry = dict(source=item["path"], unpacked=str(dest.relative_to(DATA)), bytes=dest.stat().st_size)
    with h5py.File(dest, "r") as f:
        matrix = f["0/DATA/0/matrix"]
        entry.update(matrix_shape=list(matrix.shape), dtype=str(matrix.dtype),
                     gene_id_count=len(f["0/META/ROW/id"]), signature_id_count=len(f["0/META/COL/id"]))
        assert matrix.size == entry["gene_id_count"] * entry["signature_id_count"]
        matrix[0, :4]
    records.append(entry)
    print(json.dumps(entry, ensure_ascii=False), flush=True)
out = DATA / "unpacked/LINCS_matrix_contents.json"
out.write_text(json.dumps(records, ensure_ascii=False, indent=2))
