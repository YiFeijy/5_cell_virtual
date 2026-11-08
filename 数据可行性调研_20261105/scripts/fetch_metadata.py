#!/usr/bin/env python3
"""Re-fetch the audited bounded source set with curl; no credentials required.

Existing validated files are retained. A .part download is only committed after
format validation, preventing a login page from becoming a CSV/ZIP/Parquet file.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def valid(path, expected):
    if not path.exists() or not path.stat().st_size:
        return False
    head = path.read_bytes()[:64].lstrip().lower()
    if expected.endswith((".zip", ".whl", ".xlsx")):
        return zipfile.is_zipfile(path)
    if expected.endswith(".parquet"):
        return head.startswith(b"par1")
    if expected.endswith(".gz"):
        return head.startswith(b"\x1f\x8b")
    if expected.endswith(".xls"):
        return head.startswith(b"\xd0\xcf\x11\xe0")
    if expected.endswith((".csv", ".tsv", ".rds", ".ipynb", ".json", ".md")):
        return not head.startswith((b"<!doctype", b"<html"))
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-expression", action="store_true", help="Also fetch the four original MIX-seq expression ZIPs (~103 MB total)")
    args = parser.parse_args()
    rows = list(csv.DictReader((ROOT / "download_sources.csv").open()))
    results = []
    for row in rows:
        dest = ROOT / row["path"]
        if "mixseq_" in dest.name and dest.suffix == ".zip" and not args.include_expression:
            continue
        if valid(dest, row["path"]):
            status = "existing_format_valid"
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            part = dest.with_name(dest.name + ".part")
            proc = subprocess.run(["curl", "-fsSL", "--max-time", "90", "--max-filesize", "125829120", "-o", str(part), row["url"]], capture_output=True, text=True)
            if proc.returncode or not valid(part, row["path"]):
                results.append(dict(row, status="failed_or_nondata_response", details=proc.stderr.strip()[:300]))
                if part.exists():
                    part.unlink()
                continue
            part.replace(dest)
            status = "downloaded_format_valid"
        results.append(dict(row, status=status, bytes=dest.stat().st_size, sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
    keys = ["path", "url", "status", "bytes", "sha256", "details"]
    with (ROOT / "fetch_log.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)
    print("Validated", sum(x["status"] != "failed_or_nondata_response" for x in results), "of", len(results), "sources; see fetch_log.csv")


if __name__ == "__main__":
    main()
