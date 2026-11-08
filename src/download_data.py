#!/usr/bin/env python3
"""Resumable public-data downloads, official checksums and a live purpose table."""
import argparse
import concurrent.futures
import csv
from datetime import datetime
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import threading
import time
import zipfile
from zoneinfo import ZoneInfo

DATA = Path(__file__).resolve().parents[1] / 'data'
LOCK = threading.RLock()
STATES = {}
ROWS = []
PREVIOUS = {}


def now():
    return datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="seconds")


def update(row, **values):
    with LOCK:
        STATES[row["id"]].update(values, updated_at=now())


def digest(path, algorithm):
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        while data := f.read(8 * 1024 * 1024):
            h.update(data)
    return h.hexdigest()


def verify(row, path):
    size = path.stat().st_size
    expected = row.get("expected_size")
    if expected and size != expected:
        raise ValueError(f"size {size} != official size {expected}")
    with path.open("rb") as f:
        head = f.read(128).lstrip().lower()
    if head.startswith((b"<html", b"<!doctype")):
        raise ValueError("received HTML instead of a data file")
    checksum = row.get("checksum")
    algorithm, wanted = checksum.split(":", 1) if checksum else ("sha256", None)
    got = digest(path, algorithm)
    if wanted and got.lower() != wanted.lower():
        raise ValueError(f"official {algorithm} mismatch")
    name = row["name"]
    details = {"bytes": size, "checksum_algorithm": algorithm, "checksum": got,
               "official_checksum_verified": bool(wanted), "verified_mtime_ns": path.stat().st_mtime_ns}
    if name.endswith((".zip", ".xlsx")):
        with zipfile.ZipFile(path) as z:
            broken = z.testzip()
            if broken:
                raise ValueError(f"ZIP CRC failed: {broken}")
            details["archive_members"] = len(z.infolist())
    elif name.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            details["decompressed_header_hex"] = f.read(32).hex()
        # Published SHA512 verifies complete compressed matrices; no 50GB expansion.
        if not wanted:
            with gzip.open(path, "rb") as f:
                while f.read(8 * 1024 * 1024):
                    pass
    elif name.endswith((".parquet", ".h5ad")):
        if name.endswith(".parquet"):
            import pyarrow.parquet as pq
            m = pq.read_metadata(path)
            details.update(rows=m.num_rows, columns=m.num_columns)
        else:
            import h5py
            with h5py.File(path, "r") as f:
                assert "X" in f and "obs" in f and "var" in f, list(f)
                shape = f["X"].attrs.get("shape")
                if shape is None and isinstance(f["X"], h5py.Dataset):
                    shape = f["X"].shape
                if shape is not None:
                    details["expression_shape"] = [int(x) for x in shape]
    elif name.endswith(".csv"):
        with path.open(encoding="utf-8-sig", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
            next(reader)
        details["csv_columns"] = len(header)
    elif name.endswith(".xls"):
        if not head.startswith(b"\xd0\xcf\x11\xe0"):
            raise ValueError("not an XLS workbook")
    return details


def worker(row):
    path = DATA / row["path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    logpath = DATA / "download_logs" / (row["id"].replace("/", "__") + ".log")
    try:
        previous = PREVIOUS.get(row["id"], {})
        previous_digest = f"{previous.get('checksum_algorithm')}:{previous.get('checksum')}"
        if (path.exists() and previous.get("status") == "已完成"
                and path.stat().st_size == previous.get("bytes")
                and (not row.get("checksum") or row["checksum"] == previous_digest)
                and (not previous.get("verified_mtime_ns") or path.stat().st_mtime_ns == previous["verified_mtime_ns"])):
            info = {key: previous[key] for key in ["bytes", "checksum_algorithm", "checksum", "official_checksum_verified", "verified_mtime_ns", "archive_members", "decompressed_header_hex", "rows", "columns", "expression_shape", "csv_columns"] if key in previous}
            update(row, status="已完成", verification_reused=True, **info)
            return
        if path.exists():
            update(row, status="校验中")
            info = verify(row, path)
            update(row, status="已完成", **info)
            return
        # A prior transfer may have finished before the manager took it over.
        # Resuming at EOF can return HTTP 416 even when the file is complete.
        if part.exists() and ((row.get("expected_size") and part.stat().st_size == row["expected_size"])
                              or (row["name"].endswith((".zip", ".xlsx")) and zipfile.is_zipfile(part))):
            update(row, status="校验中")
            info = verify(row, part)
            part.replace(path)
            update(row, status="已完成", **info)
            return
        update(row, status="下载中", started_at=now())
        # Start a fresh curl process after each interrupted connection so -C -
        # always uses the current file length. A time-limited transfer that made
        # progress is resumed, rather than dropping a multi-GB file permanently.
        for attempt in range(1, 97):
            before = part.stat().st_size if part.exists() else 0
            update(row, status="下载中", attempt=attempt)
            with logpath.open("a") as log:
                command = ["curl", "-fsSL", "--retry", "0", "--connect-timeout", "30", "--max-time", "600",
                           "--speed-time", "120", "--speed-limit", "1024", "-C", "-", "-o", str(part), row["url"]]
                log.write(f"\n[{now()}] Attempt {attempt}, resume at {before}: {row['url']}\n"); log.flush()
                proc = subprocess.run(command, stdout=log, stderr=log)
            after = part.stat().st_size if part.exists() else 0
            if not proc.returncode or (row.get("expected_size") and after == row["expected_size"]):
                break
            if proc.returncode == 33:
                raise RuntimeError("Server refused range resume; partial file retained")
            if attempt == 96:
                raise RuntimeError(f"curl exit {proc.returncode}; retry limit reached; .part retained")
            update(row, status="等待续传", last_curl_exit=proc.returncode,
                   last_attempt_progress_bytes=after-before)
            time.sleep(min(60, attempt * 5))
        update(row, status="校验中")
        info = verify(row, part)
        part.replace(path)
        update(row, status="已完成", **info)
        print(now(), "COMPLETE", row["id"], info["bytes"], flush=True)
    except Exception as e:
        update(row, status="失败", error=str(e))
        print(now(), "FAILED", row["id"], str(e), flush=True)


def render():
    with LOCK:
        snapshot = [dict(STATES[row["id"]]) for row in ROWS]
    for item in snapshot:
        path = DATA / item["path"]
        part = path.with_name(path.name + ".part")
        item["downloaded_bytes"] = path.stat().st_size if path.exists() else (part.stat().st_size if part.exists() else 0)
    obj = {"updated_at": now(), "files": snapshot}
    tmp = DATA / "download_status.json.tmp"
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2))
    tmp.replace(DATA / "download_status.json")
    cols = ["dataset", "name", "role", "purpose", "status", "downloaded_bytes", "expected_size", "path", "scope", "error"]
    with (DATA / "download_status.csv.tmp").open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(snapshot)
    (DATA / "download_status.csv.tmp").replace(DATA / "download_status.csv")
    complete = [x for x in snapshot if x["status"] == "已完成"]
    text = ["# 数据用途与下载状态", "", f"更新时间：{obj['updated_at']}。已完成 {len(complete)}/{len(snapshot)} 个文件，完整文件共 {sum(x['downloaded_bytes'] for x in complete)/1e9:.2f} GB。", "",
            "先保存数据并检查文件，再讨论对接和模型。以下逐文件写明内容用途；`.part`为未完成文件。", "",
            "| 数据集 | 文件 | 数据是干什么的 | 状态 | 已落盘大小 |", "|---|---|---|---|---|"]
    for item in snapshot:
        display = f"{item['downloaded_bytes']/1e6:.1f} MB"
        if item.get("expected_size"):
            display += f" / {item['expected_size']/1e6:.1f} MB"
        text.append(f"| {item['dataset']} | [{item['name']}]({item['path']}) | {item['purpose']} | {item['status']} | {display} |")
    text += ["", "说明：CTRP PharmacoSet是整理版，原始发布ZIP另存；药敏曲线汇总和浓度级实测分别保存。Tahoe本批只有一个完整表达分片及三张说明表，不代表全库下载完成。", "",
             "文件完成后检查官方校验码（来源有提供时）、文件大小及格式；ZIP检查CRC。SHA512已核验的LINCS大矩阵暂不展开，以保留原下载文件。", ""]
    tmp = DATA / "数据用途与下载状态.md.tmp"
    tmp.write_text("\n".join(text)); tmp.replace(DATA / "数据用途与下载状态.md")
    return snapshot


def main():
    global ROWS, PREVIOUS
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--max-priority", type=int, default=100)
    parser.add_argument("--dataset", action="append", help="Only download this dataset; repeat to select several datasets.")
    parser.add_argument("--dry-run", action="store_true", help="List selected files without downloading or changing status.")
    args = parser.parse_args()
    ROWS = [r for r in json.loads((DATA / "download_manifest.json").read_text()) if r["priority"] <= args.max_priority]
    if args.dataset:
        available = {r['dataset'] for r in ROWS}
        unknown = set(args.dataset) - available
        if unknown:
            parser.error(f"Unknown or priority-excluded dataset(s): {sorted(unknown)}; available: {sorted(available)}")
        ROWS = [r for r in ROWS if r['dataset'] in set(args.dataset)]
    if args.dry_run:
        print(json.dumps({'selected_files': len(ROWS), 'files': ROWS}, ensure_ascii=False, indent=2))
        return
    (DATA / 'download_logs').mkdir(parents=True, exist_ok=True)
    lock = (DATA / "download.lock").open("w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Download manager is already running")
    statusfile = DATA / "download_status.json"
    if statusfile.exists():
        PREVIOUS = {x["id"]: x for x in json.loads(statusfile.read_text()).get("files", [])}
    for row in ROWS:
        STATES[row["id"]] = dict(row, status="排队中", updated_at=now())
    render()
    print(now(), "START", len(ROWS), "files", "free_GB", round(shutil.disk_usage(DATA).free/1e9, 1), flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(worker, row) for row in ROWS]
        while not all(f.done() for f in futures):
            render()
            time.sleep(5)
        for f in futures:
            f.result()
    snapshot = render()
    failed = [x["id"] for x in snapshot if x["status"] == "失败"]
    print(now(), "FINISHED", "failed", failed, flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
