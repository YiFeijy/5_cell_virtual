#!/usr/bin/env python3
"""Unpack complete functional datasets and keep original archives unchanged."""
from pathlib import Path
import json
import shutil
import zipfile

DATA = Path(__file__).resolve().parents[1]
sources = [
    ("raw/CTRPv2/CTRPv2.0_2015_ctd2_ExpandedDataset.zip", "CTRPv2_2015"),
    ("raw/CTRPv2/CTRPv2.1_2016_pub_NatChemBiol_12_109.zip", "CTRPv2_2016"),
    ("raw/GDSC/GDSC1_public_raw_data_27Oct23.zip", "GDSC1"),
    ("raw/GDSC/GDSC2_public_raw_data_27Oct23.zip", "GDSC2"),
]
records = []
for relative, group in sources:
    source = DATA / relative
    dest = DATA / "unpacked" / group
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source) as z:
        for member in z.infolist():
            target = dest / member.filename
            if not target.resolve().is_relative_to(dest.resolve()):
                raise ValueError(f"Unsafe archive member: {member.filename}")
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or target.stat().st_size != member.file_size:
                part = target.with_name(target.name + ".unpacking")
                with z.open(member) as reader, part.open("wb") as writer:
                    shutil.copyfileobj(reader, writer, length=8 * 1024 * 1024)
                assert part.stat().st_size == member.file_size
                part.replace(target)
            records.append(dict(source=relative, file=str(target.relative_to(DATA)), bytes=member.file_size))
    print("Unpacked", group, flush=True)
(DATA / "unpacked/manifest.json").write_text(json.dumps(records, ensure_ascii=False, indent=2))
print("Unpacked files:", len(records), "; datasets were not joined")
