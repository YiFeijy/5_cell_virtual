#!/usr/bin/env python3
"""
Step 00: 原始下载数据体检与结构剖析脚本（可手动复现）

用途：
1. 统计本地已下载的全部原始数据集（文件清单、格式、大小）；
2. 深入检查最初版本 (v0) 必需的三大核心原始数据源：
   - HMS_LINCS (Screen20237/20238/20239 原孔 Excel)
   - LINCS_GSE70138 (meta 表 + Level3/Level5 GCTX 矩阵)
   - CCLE_DepMap_24Q4 (Model.csv + 19,193 基因基线表达矩阵)
3. 概览检查后续迭代将用到的扩展数据集（LINCS_2020、CTRPv2、CellModelPassports 等）；
4. 将实测结果落盘为 raw_data_inventory.csv 和 raw_data_profile.json，供后续预处理和建模核查。

用法示例：
  python iterations/00_raw_data_audit/inspect_raw_data.py
  python iterations/00_raw_data_audit/inspect_raw_data.py --raw-root /path/to/data/raw --unpacked-root /path/to/data/unpacked
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent


def resolve_default_roots() -> tuple[Path, Path | None]:
    """优先使用上级工作区 3_MVP/data/raw（含全部 134 GB 原始大文件），否则退回仓库内 data/raw。"""
    workspace_raw = REPO_ROOT.parent / "data" / "raw"
    workspace_unpacked = REPO_ROOT.parent / "data" / "unpacked"
    repo_raw = REPO_ROOT / "data" / "raw"
    repo_unpacked = REPO_ROOT / "data" / "unpacked"
    raw_root = workspace_raw if workspace_raw.exists() else repo_raw
    if workspace_unpacked.exists():
        unpacked_root = workspace_unpacked
    elif repo_unpacked.exists():
        unpacked_root = repo_unpacked
    else:
        unpacked_root = None
    return raw_root, unpacked_root


def human_size(n_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    v = float(n_bytes)
    for u in units:
        if v < 1024.0 or u == units[-1]:
            return f"{v:.2f} {u}"
        v /= 1024.0
    return f"{n_bytes} B"


def inventory_raw_dir(raw_root: Path) -> tuple[list[dict], list[dict]]:
    """逐目录、逐文件统计已下载原始数据。"""
    file_rows = []
    dir_summary = []
    for sub in sorted(p for p in raw_root.iterdir() if p.is_dir()):
        files = sorted(f for f in sub.rglob("*") if f.is_file())
        total_bytes = sum(f.stat().st_size for f in files)
        exts = sorted({"".join(f.suffixes[-2:]) if f.name.endswith(".gz") else f.suffix for f in files})
        dir_summary.append({
            "dataset": sub.name,
            "n_files": len(files),
            "total_bytes": total_bytes,
            "total_size_human": human_size(total_bytes),
            "file_types": ", ".join(exts),
        })
        for f in files:
            file_rows.append({
                "dataset": sub.name,
                "relative_path": str(f.relative_to(raw_root)),
                "bytes": f.stat().st_size,
                "size_human": human_size(f.stat().st_size),
            })
    return dir_summary, file_rows


def inspect_hms(raw_root: Path) -> dict:
    """深入剖析 HMS LINCS 三批原孔 Excel（最初版本功能标签的唯一来源）。"""
    hms_dir = raw_root / "HMS_LINCS"
    if not hms_dir.exists():
        hms_dir = REPO_ROOT / "data" / "raw" / "HMS_LINCS"
    xls_files = sorted(hms_dir.glob("*.xls"))
    if not xls_files:
        return {"status": "missing"}

    frames = []
    per_file = []
    sheets_example = None
    for path in xls_files:
        xl = pd.ExcelFile(path, engine="xlrd")
        if sheets_example is None:
            sheets_example = xl.sheet_names
        df = xl.parse("Data")
        df["source_file"] = path.name
        frames.append(df)
        per_file.append({
            "file": path.name,
            "bytes": path.stat().st_size,
            "sheets": xl.sheet_names,
            "n_wells": int(len(df)),
        })

    raw = pd.concat(frames, ignore_index=True)
    rcc = pd.to_numeric(raw["Relative Cell Count"], errors="coerce")
    rgr = pd.to_numeric(raw["Relative Growth Rate After Treatment"], errors="coerce")
    n0 = pd.to_numeric(raw["Total Cell Count Before Treatment"], errors="coerce")
    nt = pd.to_numeric(raw["Total Cell Count After Treatment"], errors="coerce")
    nc = pd.to_numeric(raw["Total Control Cell Count"], errors="coerce")
    doses = sorted(pd.to_numeric(raw["Small Mol Concentration"], errors="coerce").dropna().unique().tolist())
    cells = sorted(raw["CellLine"].dropna().astype(str).unique().tolist())
    plates = sorted(raw["LJP Library Plate"].dropna().astype(str).unique().tolist())

    # 抽取经典示范样本：BT20 + afatinib (BRD-A58767537) + 0.04 uM
    ex_mask = (
        (raw["CellLine"].astype(str) == "BT20")
        & (raw["BRDid"].astype(str) == "BRD-A58767537")
        & np.isclose(pd.to_numeric(raw["Small Mol Concentration"], errors="coerce"), 0.04)
    )
    ex_cols = [
        "source_file", "LJP Library Plate", "CellLine", "BRDid",
        "Small Mol Concentration", "Timepoint",
        "Total Cell Count Before Treatment", "Total Cell Count After Treatment",
        "Total Control Cell Count", "Relative Cell Count",
        "Relative Growth Rate After Treatment",
    ]
    examples = raw.loc[ex_mask, [c for c in ex_cols if c in raw.columns]].head(3).to_dict(orient="records")

    return {
        "status": "ok",
        "files": per_file,
        "total_wells": int(len(raw)),
        "columns": raw.columns.tolist(),
        "cells": cells,
        "n_cells": len(cells),
        "n_unique_drugs_BRDid": int(raw["BRDid"].nunique()),
        "doses_uM": doses,
        "timepoints_h": sorted(raw["Timepoint"].dropna().unique().tolist()),
        "conc_units": sorted(raw["Small Mol Conc Unit"].dropna().astype(str).unique().tolist()),
        "library_plates": plates,
        "value_stats": {
            "N0_before_treatment": {"min": float(n0.min()), "median": float(n0.median()), "max": float(n0.max())},
            "Nt_after_treatment": {"min": float(nt.min()), "median": float(nt.median()), "max": float(nt.max())},
            "Nc_control_72h": {"min": float(nc.min()), "median": float(nc.median()), "max": float(nc.max())},
            "relative_cell_count_Nt_over_Nc": {
                "min": round(float(rcc.min()), 4),
                "q25": round(float(rcc.quantile(0.25)), 4),
                "median": round(float(rcc.median()), 4),
                "q75": round(float(rcc.quantile(0.75)), 4),
                "max": round(float(rcc.max()), 4),
            },
            "relative_net_growth_Nt_minus_N0_over_Nc_minus_N0": {
                "min": round(float(rgr.min()), 4),
                "q25": round(float(rgr.quantile(0.25)), 4),
                "median": round(float(rgr.median()), 4),
                "q75": round(float(rgr.quantile(0.75)), 4),
                "max": round(float(rgr.max()), 4),
                "negative_well_count_net_killing": int((rgr < 0).sum()),
                "negative_well_fraction": round(float((rgr < 0).mean()), 4),
            },
        },
        "example_wells_BT20_afatinib_0.04uM": examples,
    }


def inspect_lincs_gse70138(raw_root: Path, unpacked_root: Path | None) -> dict:
    """深入剖析 LINCS Phase II (GSE70138) 元数据与 GCTX 表达矩阵。"""
    lincs_dir = raw_root / "LINCS_GSE70138"
    if not lincs_dir.exists():
        lincs_dir = REPO_ROOT / "data" / "raw" / "LINCS_GSE70138"
    if not lincs_dir.exists():
        return {"status": "missing"}

    gene_df = pd.read_csv(next(lincs_dir.glob("*gene_info*")), sep="\t")
    cell_df = pd.read_csv(next(lincs_dir.glob("*cell_info*")), sep="\t")
    pert_df = pd.read_csv(next(lincs_dir.glob("*pert_info*")), sep="\t")
    sig_df = pd.read_csv(next(lincs_dir.glob("*sig_info*")), sep="\t", low_memory=False)
    inst_df = pd.read_csv(next(lincs_dir.glob("*inst_info*")), sep="\t", low_memory=False)

    trt_cp_sig = sig_df[sig_df["pert_type"] == "trt_cp"]
    ljp_sig = trt_cp_sig[trt_cp_sig["sig_id"].str.startswith(("LJP005_", "LJP006_"))]

    # 统计 canonical_smiles 的缺失或受限情况（例如 restricted / -666）
    pert_cp = pert_df[pert_df["pert_type"] == "trt_cp"]
    smi_s = pert_cp["canonical_smiles"].fillna("").astype(str).str.strip()
    smi_restricted = int(smi_s.str.lower().isin({"restricted", "-666", ""}).sum())

    res = {
        "status": "ok",
        "gene_info": {
            "total_genes": int(len(gene_df)),
            "measured_landmark_pr_is_lm_1": int(gene_df["pr_is_lm"].eq(1).sum()),
            "inferred_genes": int(gene_df["pr_is_lm"].eq(0).sum()),
            "columns": gene_df.columns.tolist(),
        },
        "cell_info": {
            "rows_in_cell_info_table": int(len(cell_df)),
            "unique_cell_ids_observed_in_sig_info": int(sig_df["cell_id"].nunique()),
            "unique_cell_ids_with_trt_cp": int(trt_cp_sig["cell_id"].nunique()),
            "hms_6_cells_present_in_trt_cp": sorted(
                set(trt_cp_sig["cell_id"]) & {"BT20", "HS578T", "MCF10A", "MCF7", "MDAMB231", "SKBR3"}
            ),
        },
        "pert_info": {
            "total_perts": int(len(pert_df)),
            "trt_cp_perts": int(len(pert_cp)),
            "trt_cp_restricted_or_missing_smiles": smi_restricted,
            "columns": pert_df.columns.tolist(),
        },
        "sig_info": {
            "total_signatures": int(len(sig_df)),
            "pert_type_counts": sig_df["pert_type"].value_counts().to_dict(),
            "trt_cp_signatures": int(len(trt_cp_sig)),
            "trt_cp_time_distribution": trt_cp_sig["pert_itime"].value_counts().head(10).to_dict(),
            "ljp005_006_trt_cp_signatures": int(len(ljp_sig)),
            "ljp005_006_time_distribution": ljp_sig["pert_itime"].value_counts().to_dict(),
            "columns": sig_df.columns.tolist(),
        },
        "inst_info": {
            "total_instances": int(len(inst_df)),
            "pert_type_counts": inst_df["pert_type"].value_counts().head(8).to_dict(),
            "ctl_vehicle_instances": int((inst_df["pert_type"] == "ctl_vehicle").sum()),
            "columns": inst_df.columns.tolist(),
        },
        "unpacked_gctx": {},
    }

    if unpacked_root is not None and (unpacked_root / "LINCS_GSE70138").exists():
        for tag, pat in [("Level3", "*Level3*.gctx"), ("Level5", "*Level5*.gctx")]:
            matches = list((unpacked_root / "LINCS_GSE70138").glob(pat))
            if matches:
                p = matches[0]
                with h5py.File(p, "r") as f:
                    mat = f["0/DATA/0/matrix"]
                    sample_block = mat[:64, :978]
                    res["unpacked_gctx"][tag] = {
                        "file": p.name,
                        "bytes": p.stat().st_size,
                        "size_human": human_size(p.stat().st_size),
                        "shape_rows_x_genes": [int(x) for x in mat.shape],
                        "dtype": str(mat.dtype),
                        "sample_64x978_min": round(float(np.min(sample_block)), 4),
                        "sample_64x978_median": round(float(np.median(sample_block)), 4),
                        "sample_64x978_max": round(float(np.max(sample_block)), 4),
                    }
    return res


def inspect_ccle(raw_root: Path) -> dict:
    """深入剖析 CCLE / DepMap 24Q4 细胞模型表与基线表达矩阵。"""
    ccle_dir = raw_root / "CCLE_DepMap_24Q4"
    if not ccle_dir.exists():
        return {"status": "missing"}

    model_df = pd.read_csv(ccle_dir / "Model.csv")
    expr_path = ccle_dir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"

    # 轻量读取表头 + 第一列 ModelID（避免不必要的全矩阵内存开销）
    with expr_path.open("r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        expr_model_ids = []
        sample_vals = []
        for idx, row in enumerate(reader):
            if not row:
                continue
            expr_model_ids.append(row[0])
            if idx < 5:
                sample_vals.extend(float(v) for v in row[1:201])
    expr_id_set = set(expr_model_ids)

    # 核查 HMS 6 株细胞在 CCLE 24Q4 的身份与 RNA 覆盖（解释为什么 MCF10A 进不了 v0）
    hms_targets = ["BT20", "HS578T", "MCF10A", "MCF7", "MDAMB231", "SKBR3"]
    stripped_map = model_df.assign(
        _key=model_df["StrippedCellLineName"].fillna("").astype(str).str.upper()
    ).set_index("_key")
    hms_check = []
    for cname in hms_targets:
        if cname in stripped_map.index:
            rec = stripped_map.loc[cname]
            if isinstance(rec, pd.DataFrame):
                rec = rec.iloc[0]
            mid = str(rec["ModelID"])
            hms_check.append({
                "hms_cell_name": cname,
                "ModelID": mid,
                "OncotreeLineage": str(rec.get("OncotreeLineage", "")),
                "OncotreePrimaryDisease": str(rec.get("OncotreePrimaryDisease", "")),
                "has_24Q4_RNA_row": mid in expr_id_set,
            })
        else:
            hms_check.append({"hms_cell_name": cname, "ModelID": None, "has_24Q4_RNA_row": False})

    return {
        "status": "ok",
        "model_csv_rows": int(len(model_df)),
        "model_csv_key_columns": ["ModelID", "CellLineName", "StrippedCellLineName", "RRID", "OncotreeLineage", "SangerModelID", "COSMICID"],
        "expression_matrix": {
            "file": expr_path.name,
            "bytes": expr_path.stat().st_size,
            "size_human": human_size(expr_path.stat().st_size),
            "n_models_with_RNA": len(expr_model_ids),
            "n_protein_coding_genes": len(header) - 1,
            "example_gene_columns": header[1:6],
            "unit": "log2(TPM + 1)",
            "sample_5x200_min": round(float(min(sample_vals)), 4),
            "sample_5x200_median": round(float(np.median(sample_vals)), 4),
            "sample_5x200_max": round(float(max(sample_vals)), 4),
        },
        "hms_6_cells_coverage_audit": hms_check,
    }


def inspect_extension_datasets(raw_root: Path, unpacked_root: Path | None) -> dict:
    """对后续迭代（v1/v2/v3 及外部盲测）将使用的扩展数据做快速实测概览。"""
    ext = {}

    # 1. LINCS_2020 (v1/v2/v3 主力转录扩展)
    l20 = raw_root / "LINCS_2020"
    if (l20 / "siginfo_beta.txt").exists():
        sig20 = pd.read_csv(l20 / "siginfo_beta.txt", sep="\t", usecols=["pert_type", "cell_iname", "pert_time", "qc_pass"], low_memory=False)
        trt = sig20[sig20["pert_type"] == "trt_cp"]
        comp20 = pd.read_csv(l20 / "compoundinfo_beta.txt", sep="\t", usecols=["pert_id", "canonical_smiles"], low_memory=False)
        ext["LINCS_2020"] = {
            "role_in_roadmap": "v1–v3 大规模小分子多时间点转录监督（核心化学泛化杠杆）",
            "total_trt_cp_signatures": int(len(trt)),
            "qc_pass_trt_cp_signatures": int((trt["qc_pass"] == 1).sum()),
            "unique_cells_trt_cp": int(trt["cell_iname"].nunique()),
            "unique_pert_ids_in_compoundinfo": int(comp20["pert_id"].nunique()),
            "top_timepoints_h": trt["pert_time"].value_counts().head(6).to_dict(),
        }

    # 2. CTRPv2 (v1/v2/v3 多谱系 72h 药敏扩展)
    if unpacked_root is not None and (unpacked_root / "CTRPv2_2015" / "v20.meta.per_compound.txt").exists():
        c_cpd = pd.read_csv(unpacked_root / "CTRPv2_2015" / "v20.meta.per_compound.txt", sep="\t")
        c_ccl = pd.read_csv(unpacked_root / "CTRPv2_2015" / "v20.meta.per_cell_line.txt", sep="\t")
        ext["CTRPv2_2015"] = {
            "role_in_roadmap": "v1–v3 多谱系 72h 实测 ATP 药敏训练与筛选排序评测（剔除 49 组合，留 496 单药）",
            "compound_meta_rows": int(len(c_cpd)),
            "cell_line_meta_rows": int(len(c_ccl)),
            "dose_level_file_bytes": (unpacked_root / "CTRPv2_2015" / "v20.data.per_cpd_post_qc.txt").stat().st_size,
        }

    # 3. CellModelPassports (v3 生长率校准辅助监督)
    cmp_gr = raw_root / "CellModelPassports" / "growth_rate_20220907.csv"
    if cmp_gr.exists():
        gr = pd.read_csv(cmp_gr)
        ext["CellModelPassports"] = {
            "role_in_roadmap": "v3 训练期细胞固有 72h 对照倍增倍数辅助校准（绝不作测试必填输入）",
            "rows": int(len(gr)),
            "unique_models": int(gr["model_id"].nunique()),
            "columns": gr.columns.tolist(),
        }

    # 4. 其余已下载外部验证集清单
    for name, role in [
        ("GDSC", "外部 72h 跨测定药敏验证（GDSC1 染料 / GDSC2 ATP）"),
        ("PRISM_19Q4", "外部 120h 混合培养条形码药敏盲测（~900 非 CTRP 跨模态新药）"),
        ("PRISM_24Q2", "外部 120h 新版 PRISM 复验与 72h→120h 标度检验"),
        ("NCI_TPW", "外部跨平台（微阵列）60 细胞 × 15 药 × 2/6/24h 时间动态盲测"),
        ("Novartis_DRUGseq", "外部 U2OS × 4,343 药跨平台 MoA 检索验证（不含实测活性标签）"),
        ("MIXseq", "单细胞药物响应异质性与 3–48h trametinib 时序验证"),
        ("sciPlex", "单细胞核转录组响应与群体组成解耦（留工作二）"),
        ("Tahoe_100M", "单细胞大规模扰动分片与元数据（留工作二）"),
    ]:
        d = raw_root / name
        if d.exists():
            files = [f for f in d.rglob("*") if f.is_file()]
            ext[name] = {
                "role_in_roadmap": role,
                "n_files": len(files),
                "total_size_human": human_size(sum(f.stat().st_size for f in files)),
            }
    return ext


def main():
    ap = argparse.ArgumentParser(description="Inspect all downloaded raw datasets and profile core v0 inputs.")
    default_raw, default_unpacked = resolve_default_roots()
    ap.add_argument("--raw-root", type=Path, default=default_raw)
    ap.add_argument("--unpacked-root", type=Path, default=default_unpacked)
    args = ap.parse_args()

    print(f"=== Step 00: 原始数据体检 ===")
    print(f"raw_root      = {args.raw_root}")
    print(f"unpacked_root = {args.unpacked_root}")

    dir_summary, file_rows = inventory_raw_dir(args.raw_root)
    total_bytes = sum(d["total_bytes"] for d in dir_summary)
    print(f"\n[1/4] 本地原始数据目录数: {len(dir_summary)}, 总文件数: {len(file_rows)}, 总大小: {human_size(total_bytes)}")
    for d in dir_summary:
        print(f"  - {d['dataset']:22s} | {d['n_files']:3d} 文件 | {d['total_size_human']:>10s} | {d['file_types']}")

    print("\n[2/4] 剖析核心 v0 数据源: HMS_LINCS / LINCS_GSE70138 / CCLE_DepMap_24Q4 ...")
    hms_info = inspect_hms(args.raw_root)
    lincs_info = inspect_lincs_gse70138(args.raw_root, args.unpacked_root)
    ccle_info = inspect_ccle(args.raw_root)

    print("\n[3/4] 概览后续迭代扩展数据集 ...")
    ext_info = inspect_extension_datasets(args.raw_root, args.unpacked_root)

    profile = {
        "raw_root_scanned": str(args.raw_root),
        "total_datasets": len(dir_summary),
        "total_files": len(file_rows),
        "total_bytes": total_bytes,
        "total_size_human": human_size(total_bytes),
        "directory_summary": dir_summary,
        "v0_core_datasets": {
            "HMS_LINCS": hms_info,
            "LINCS_GSE70138": lincs_info,
            "CCLE_DepMap_24Q4": ccle_info,
        },
        "extension_datasets_for_later_iterations": ext_info,
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dir_summary).to_csv(OUT_DIR / "raw_data_inventory.csv", index=False)
    (OUT_DIR / "raw_data_profile.json").write_text(json.dumps(profile, ensure_ascii=False, indent=2) + "\n")

    print(f"\n[4/4] 完成！已保存:")
    print(f"  - {OUT_DIR / 'raw_data_inventory.csv'}")
    print(f"  - {OUT_DIR / 'raw_data_profile.json'}")
    if hms_info.get("status") == "ok":
        print(f"\n[核心核对摘要]")
        print(f"  HMS_LINCS: {hms_info['total_wells']} 原孔, {hms_info['n_cells']} 细胞 {hms_info['cells']}, {hms_info['n_unique_drugs_BRDid']} 药, 净杀伤(rgr<0)占比 {hms_info['value_stats']['relative_net_growth_Nt_minus_N0_over_Nc_minus_N0']['negative_well_fraction']*100:.1f}%")
    if lincs_info.get("status") == "ok":
        print(f"  LINCS_GSE70138: {lincs_info['sig_info']['total_signatures']} 签名 ({lincs_info['sig_info']['trt_cp_signatures']} trt_cp), {lincs_info['gene_info']['measured_landmark_pr_is_lm_1']}/{lincs_info['gene_info']['total_genes']} 实测 landmark 基因")
    if ccle_info.get("status") == "ok":
        print(f"  CCLE_DepMap_24Q4: {ccle_info['expression_matrix']['n_models_with_RNA']} 个 RNA 细胞模型 × {ccle_info['expression_matrix']['n_protein_coding_genes']} 基因; HMS 细胞命中:")
        for c in ccle_info["hms_6_cells_coverage_audit"]:
            print(f"    * {c['hms_cell_name']:8s} -> {c['ModelID']} | has_24Q4_RNA={c['has_24Q4_RNA_row']}")


if __name__ == "__main__":
    main()
