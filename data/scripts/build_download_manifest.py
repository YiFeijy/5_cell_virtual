#!/usr/bin/env python3
"""Create a source-backed download queue and copy already downloaded inputs."""
import csv
import gzip
import json
from pathlib import Path
import shutil
import re
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OLD = ROOT / "数据可行性调研_20261105"
rows = []


def add(dataset, name, url, role, purpose, size=None, checksum=None, priority=10, local=None, scope="完整文件"):
    dest = DATA / "raw" / dataset / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Adopt earlier direct downloads, including resumable partial files.
    for old in [DATA / "raw" / name, DATA / "raw" / (name + ".part")]:
        target = dest if old.name == name else dest.with_name(name + ".part")
        if old.exists() and old != target and not target.exists():
            old.rename(target)
    if local and Path(local).exists() and not dest.exists():
        shutil.copy2(local, dest)
    rows.append(dict(id=f"{dataset}/{name}", dataset=dataset, name=name, url=url,
                     path=str(dest.relative_to(DATA)), role=role, purpose=purpose,
                     expected_size=size, checksum=checksum, priority=priority, scope=scope))


def figshare_entries(filename, dataset, names, role, purpose, priority=10):
    obj = json.loads((OLD / "evidence" / filename).read_text())
    for item in obj["files"]:
        if item["name"] in names:
            md5 = item.get("computed_md5") or item.get("supplied_md5")
            add(dataset, item["name"], item["download_url"], role, purpose,
                item["size"], f"md5:{md5}" if md5 else None, priority)


figshare_entries("figshare_depmap24q4.json", "CCLE_DepMap_24Q4",
                 ["OmicsExpressionProteinCodingGenesTPMLogp1.csv"], "细胞初始状态",
                 "未给药细胞系的基线RNA表达，用于描述细胞背景", 0)
for name, old in [("Model.csv", "depmap24q4_Model.csv"),
                  ("OmicsDefaultModelProfiles.csv", "depmap24q4_OmicsDefaultModelProfiles.csv"),
                  ("OmicsProfiles.csv", "depmap24q4_OmicsProfiles.csv")]:
    obj = json.loads((OLD / "evidence/figshare_depmap24q4.json").read_text())
    item = next(x for x in obj["files"] if x["name"] == name)
    add("CCLE_DepMap_24Q4", name, item["download_url"], "细胞说明",
        "说明细胞系身份及RNA数据记录", item["size"], "md5:" + item["computed_md5"], 0,
        OLD / "metadata" / old)

figshare_entries("figshare_prism.json", "PRISM_19Q4", [
    "secondary-screen-replicate-collapsed-logfold-change.csv",
    "secondary-screen-replicate-collapsed-treatment-info.csv",
    "secondary-screen-cell-line-info.csv", "secondary-screen-readme.txt",
    "primary-screen-replicate-collapsed-logfold-change.csv",
    "primary-screen-replicate-collapsed-treatment-info.csv",
    "primary-screen-cell-line-info.csv", "primary-screen-readme.txt",
    "secondary-screen-dose-response-curve-parameters.csv"],
    "药物功能读数", "药物处理后的相对细胞丰度及药敏曲线；处理表说明药物和剂量，README说明读数定义", 1)
for old, name in [("prism_secondary_treatment.csv", "secondary-screen-replicate-collapsed-treatment-info.csv"),
                  ("prism_secondary_cell.csv", "secondary-screen-cell-line-info.csv")]:
    shutil.copy2(OLD / "metadata" / old, DATA / "raw/PRISM_19Q4" / name)
shutil.copy2(OLD / "evidence/prism_secondary_readme.txt", DATA / "raw/PRISM_19Q4/secondary-screen-readme.txt")

for version in [1, 2]:
    base = "CTRPv2.0_2015_ctd2_ExpandedDataset" if version == 1 else "CTRPv2.1_2016_pub_NatChemBiol_12_109"
    add("CTRPv2", base + ".zip", f"https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/{base}.zip",
        "药物功能读数及基线说明", "原始发布压缩包：浓度级活性、实验信息及相关基线表达资料；先保存原文件", priority=1)
c = json.loads((DATA / "manifests/ctrp_4429239.json").read_text())["files"][0]
add("CTRPv2", "CTRPv2_PharmacoSet.rds", c["links"]["self"], "药物功能读数整理版",
    "维护者整理的CTRP数据对象，包含浓度响应数组及曲线汇总；与原始发布包分别保存",
    c["size"], c["checksum"], 0)

gdsc = json.loads((DATA / "manifests/gdsc_download_files.json").read_text())
for item in gdsc["included"]:
    a = item.get("attributes", {})
    url = a.get("url") or ""
    if item["type"] == "download_file" and any(s in url for s in ["GDSC1_public_raw", "GDSC2_public_raw", "GDSC1_fitted_dose", "GDSC2_fitted_dose"]):
        raw = "public_raw" in url
        size = a.get("size")
        headers = DATA / "manifests/GDSC1_fitted_headers.txt"
        reported_size = size
        if "GDSC1_fitted_dose_response_27Oct23.xlsx" in url and headers.exists():
            total = re.search(r"content-range:\s*bytes\s+\d+-\d+/(\d+)", headers.read_text(), re.I)
            if total:
                size = int(total.group(1))
        add("GDSC", url.rsplit("/", 1)[-1], url, "药物功能读数" if raw else "药敏曲线汇总",
            "原始浓度级药物实验包，用于查看给药后功能assay读数" if raw else "拟合得到的IC50、AUC等曲线统计，和浓度级实测分别保存",
            size, priority=1)
        if size != reported_size:
            rows[-1]["source_index_reported_size"] = reported_size
            rows[-1]["size_basis"] = "Official current HTTP Content-Range; index differs by one byte"

for phase, accession in [("PhaseII", "GSE70138"), ("PhaseI", "GSE92742")]:
    checks = {}
    for line in gzip.open(DATA / f"manifests/LINCS_{accession}_SHA512SUMS.txt.gz", "rt"):
        digest, name = line.split(maxsplit=1)
        checks[name.strip().lstrip("*")] = "sha512:" + digest
    for name, digest in checks.items():
        if not any(s in name for s in ["Level3", "Level5", "inst_info", "sig_info", "sig_metrics", "gene_info", "pert_info", "cell_info"]):
            continue
        is_matrix = "Level3" in name or "Level5" in name
        priority = (70 if phase == "PhaseI" else 50) if is_matrix else 2
        role = "药后分子响应" if is_matrix else "表达实验说明"
        if "Level3" in name:
            purpose = "药物处理及对照条件的归一化基因表达，用于查看药后细胞分子状态"
        elif "Level5" in name:
            purpose = "相对对照的基因表达变化签名，描述分子响应；不含独立活性标签"
        else:
            purpose = "说明表达矩阵的细胞、药物、基因、剂量、时间、重复及质量信息"
        url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file&file={quote(name)}"
        local_names = {f"{accession}_sig_info.txt.gz": "sig_info", f"{accession}_cell_info.txt.gz": "cell_info",
                       f"{accession}_gene_info.txt.gz": "gene_info", f"{accession}_pert_info.txt.gz": "pert_info"}
        local = next((OLD / "metadata" / old for old, token in local_names.items()
                      if name == f"{accession}_Broad_LINCS_{token}" + ("_2017-03-06" if accession == "GSE70138" else "") + ".txt.gz"
                      and (OLD / "metadata" / old).exists()), None)
        size = 13540497238 if "Level3" in name and accession == "GSE70138" else (5365179698 if "Level5" in name and accession == "GSE70138" else None)
        header_kind = "L3" if "Level3" in name else ("L5" if "Level5" in name else None)
        if accession == "GSE92742" and header_kind:
            headers = DATA / "manifests" / f"LINCS_PhaseI_{header_kind}_headers.txt"
            if headers.exists():
                total = re.search(r"content-range:\s*bytes\s+\d+-\d+/(\d+)", headers.read_text(), re.I)
                if total:
                    size = int(total.group(1))
        add(f"LINCS_{accession}", name, url, role, purpose, size, digest, priority, local)

sc = json.loads((DATA / "manifests/scperturb_13350497.json").read_text())
for item in sc["files"]:
    if item["key"] not in ["SrivatsanTrapnell2020_sciplex2.h5ad", "SrivatsanTrapnell2020_sciplex3.h5ad"]:
        continue
    add("sciPlex", item["key"], item["links"]["self"], "单细胞分子响应",
        "单细胞基因表达及药物剂量条件，用于观察同一药物下细胞间响应差异", item["size"], item["checksum"], 20)

mix = json.loads((OLD / "evidence/figshare_mixseq.json").read_text())
for item in mix["files"]:
    name = item["name"]
    if name.startswith("sg"):
        continue
    md5 = item.get("computed_md5") or item.get("supplied_md5")
    old = OLD / "metadata" / ("mixseq_" + name)
    if name == "all_CL_features.rds": old = OLD / "metadata/mixseq_all_CL_features.rds"
    if name == "README.txt": old = OLD / "metadata/mixseq_README.txt"
    if name == "Supplementary Tables.xlsx": old = OLD / "metadata/mixseq_Supplementary_Tables.xlsx"
    add("MIXseq", name, item["download_url"], "单细胞分子响应及敏感性资料",
        "原始单细胞表达包、实验说明及作者提供的药敏资料；用于查看响应异质性", item["size"],
        f"md5:{md5}" if md5 else None, 10, old if old.exists() else None)

for sid, stem in [(20237, "LJP_processedData1"), (20238, "LJP_processedData2"), (20239, "LJP_processedData3")]:
    name = f"Screen{sid}_{stem}.xls"
    add("HMS_LINCS", name, f"https://lincs.hms.harvard.edu/_static/db/prod-20200624/{name}",
        "细胞生长功能读数", "处理前后细胞计数、相对细胞数和相对生长率，记录独立于表达的功能结果",
        priority=3, local=OLD / "metadata" / name)

add("NCI_TPW", "GSE116436_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116436/matrix/GSE116436_series_matrix.txt.gz",
    "时间序列表达", "NCI细胞系接受抗癌药后的多个时间点表达，用于观察分子响应随时间的变化", priority=30)

nov = json.loads((DATA / "manifests/novartis_14291446.json").read_text())
for item in nov["files"]:
    if item["key"] in ["Exp_gzip.RData", "comparisons_metadata.RData"]:
        add("Novartis_DRUGseq", item["key"], item["links"]["self"], "药后分子响应",
            "U2OS细胞接受不同药物后的表达对象及对应实验条件说明；保持官方RData格式", item["size"], item["checksum"], 35)

tahoe = json.loads((OLD / "evidence/tahoe_api.json").read_text())
one = next(x for x in tahoe if x.get("path") == "data/train-00000-of-03388.parquet")
add("Tahoe_100M", "train-00000-of-03388.parquet",
    "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/data/train-00000-of-03388.parquet",
    "单细胞分子响应", "完整的一个官方单细胞表达分片，用于确认数据内容和读取；不是Tahoe全库", one["size"], priority=15,
    scope="全库3388个表达分片中的第1个；全库未列入本批下载")
for kind in ["cell_line", "drug", "sample"]:
    add("Tahoe_100M", f"{kind}_metadata.parquet",
        f"https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/{kind}_metadata.parquet",
        "单细胞实验说明", "说明Tahoe的细胞身份、药物结构或处理样本条件", priority=2,
        local=OLD / "metadata" / f"tahoe_{kind}_metadata.parquet")

rows.sort(key=lambda x: (x["priority"], x["dataset"], x["name"]))
(DATA / "download_manifest.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2))
with (DATA / "download_manifest.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for row in rows for k in row))); w.writeheader(); w.writerows(rows)
print(f"Queued {len(rows)} files across {len(set(x['dataset'] for x in rows))} dataset folders")
