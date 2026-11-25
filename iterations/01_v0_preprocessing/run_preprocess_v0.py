#!/usr/bin/env python3
"""只读重建、核对最初 v0 队列；审计输出始终位于本迭代目录。

对应冻结的 aggregate_hms.py / extract_reference.py / prepare.py / check.py。
身份、结构分组及候选签名采用 data/reference 的固定定义，不重新定义队列。
大数组只在内存中组装；显式 --write-rebuilt-dir 才另存一份。
"""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import platform

import h5py
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DATA = REPO_ROOT / "data"
REF = DATA / "reference"
TRAIN = DATA / "training"
KEYS = ["library_id", "cell_id", "pert_id", "dose_uM"]
QC_COLUMNS = ["n_listed_instances", "n_technical_valid_instances",
              "replicate_vehicle_centered_pearson_median", "all_instances_technical_valid"]


def say(message):
    print(message, flush=True)


def records(frame):
    """让 NaN 变成严格 JSON 的 null，避免不可移植的 NaN 字面量。"""
    return json.loads(frame.to_json(orient="records", force_ascii=False, double_precision=15))


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def only_file(directory, pattern):
    paths = sorted(directory.glob(pattern))
    if len(paths) != 1:
        raise ValueError(f"需要且只允许一个 {directory / pattern}，实际找到 {len(paths)} 个")
    return paths[0]


def raw_table(kind):
    return pd.read_csv(only_file(DATA / "raw/LINCS_GSE70138", f"*{kind}*"), sep="\t")


def npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {key: z[key] for key in z.files}


def max_diff(left, right):
    a, b = np.asarray(left), np.asarray(right)
    assert a.shape == b.shape
    assert np.array_equal(np.isnan(a), np.isnan(b))
    valid = ~np.isnan(a)
    return float(np.max(np.abs(a[valid].astype(np.float64) - b[valid].astype(np.float64)), initial=0))


def exact_frame(left, right, label):
    """连同行列顺序、dtype、缺失位置及浮点数逐值完全核对。"""
    pd.testing.assert_frame_equal(left.reset_index(drop=True), right.reset_index(drop=True),
                                  check_exact=True, obj=label)
    return {"rows": len(left), "columns": len(left.columns), "exact": True,
            "numeric_max_abs_diff": 0.0}


def funnel(report, step, stage, before, frame, reason):
    row = dict(step=step, stage=stage, records_in=before, records_out=len(frame),
               records_dropped=before - len(frame),
               conditions_out=int(frame.biological_condition_group.nunique()), drop_reason=reason)
    report["funnel"].append(row)
    say(f"  {stage}: {before:,} → {len(frame):,} 条，"
        f"{row['conditions_out']:,} 个不同生物条件；{reason}")


def example_mask(frame):
    return (frame.library_id.eq("LJP006") & frame.cell_id.eq("BT20")
            & frame.pert_id.eq("BRD-A58767537") & np.isclose(frame.dose_uM, 0.04))


def stage1_hms(report):
    say("\n[Stage 1/5] HMS 原始孔 → 技术重复均值 → 三批等权均值（72 h）")
    paths = sorted((DATA / "raw/HMS_LINCS").glob("Screen2023[789]*.xls"))
    assert len(paths) == 3
    frames = []
    for path in paths:
        try:
            frame = pd.read_excel(path, sheet_name="Data", engine="xlrd")
        except ImportError as exc:
            raise SystemExit("读取 HMS .xls 需要 xlrd；请运行 pip install xlrd 后重试。") from exc
        frame["source_file"] = path.name
        frame["original_excel_row"] = np.arange(len(frame)) + 2
        frames.append(frame)
        say(f"  {path.name}: {len(frame):,} 孔")
    raw = pd.concat(frames, ignore_index=True).rename(columns={
        "CellLine": "cell_id", "BRDid": "pert_id", "Small Mol Concentration": "dose_uM",
        "Relative Cell Count": "relative_cell_count",
        "Relative Growth Rate After Treatment": "relative_net_growth"})
    assert raw.Timepoint.eq(72).all() and raw["Small Mol Conc Unit"].eq("uM").all()
    raw["library_id"] = raw["LJP Library Plate"].astype(str).str.extract(r"(\d+)")[0].astype(int).map(
        lambda value: f"LJP{value:03d}")
    raw["dose_uM"] = pd.to_numeric(raw.dose_uM)
    assert np.isfinite(raw[["relative_cell_count", "relative_net_growth"]]).all().all()
    bio = raw.groupby(KEYS + ["source_file"]).agg(
        relative_cell_count=("relative_cell_count", "mean"),
        relative_net_growth=("relative_net_growth", "mean"),
        n_technical_rows=("relative_cell_count", "size")).reset_index()
    means = bio.groupby(KEYS).agg(
        relative_cell_count=("relative_cell_count", "mean"),
        relative_net_growth=("relative_net_growth", "mean"),
        n_bio=("relative_cell_count", "size"),
        relative_cell_count_sd=("relative_cell_count", "std")).reset_index()
    means["library_condition_id"] = ["HMSv2_LIB_" + hashlib.sha256(
        f"{r.library_id}|{r.cell_id}|{r.pert_id}|{r.dose_uM:.9g}".encode()).hexdigest()[:16]
        for r in means.itertuples()]
    frozen = pd.read_csv(REF / "hms_conditions.csv")
    assert frozen.library_condition_id.is_unique and means.library_condition_id.is_unique
    matched = means.merge(frozen, on=KEYS + ["library_condition_id"], suffixes=("_raw", "_ref"),
                          validate="one_to_one")
    assert len(raw) == 36912 and len(bio) == 12420 and len(matched) == len(frozen) == 4140
    assert means.n_bio.eq(3).all() and frozen.biological_condition_group.nunique() == 3852
    diffs = {}
    for key in ["relative_cell_count", "relative_net_growth", "relative_cell_count_sd", "n_bio"]:
        diffs[key] = max_diff(matched[key + "_raw"], matched[key + "_ref"])
        assert diffs[key] < 1e-12, (key, diffs[key])
    # 生物条件不含库板：核对固定映射与实际 (cell, drug, dose) 是一一对应的。
    condition_keys = ["cell_id", "pert_id", "dose_uM"]
    assert len(frozen.drop_duplicates(condition_keys)) == 3852
    assert frozen.groupby(condition_keys).biological_condition_group.nunique().eq(1).all()
    assert frozen.groupby("biological_condition_group")[condition_keys].nunique().eq(1).all().all()
    report["stage1"] = {"source_files": [{"name": p.name, "wells": len(d), "sha256": digest(p)}
                                           for p, d in zip(paths, frames)],
                        "raw_wells": len(raw), "biorepeat_records": len(bio),
                        "library_records": len(means), "unique_conditions": 3852,
                        "time_h": 72, "dose_unit": "uM", "max_abs_diff_vs_reference": diffs,
                        "n_bio_all_three": True,
                        "technical_wells_per_biorepeat": bio.n_technical_rows.value_counts().sort_index().to_dict(),
                        "hash_rule": "HMSv2_LIB_ + sha256(library|cell|pert|dose:.9g)[:16]",
                        "assembly_label_source": "已逐条核对原孔的 hms_conditions.csv 冻结精度"}
    trace_raw = raw[example_mask(raw)]
    trace_bio = bio[example_mask(bio)]
    trace_mean = means[example_mask(means)]
    assert len(trace_mean) == 1 and len(trace_bio) == 3
    assert abs(trace_mean.relative_cell_count.iloc[0] - 0.7512411111111111) < 1e-12
    report["example_trace"] = {
        "condition": "LJP006 | BT20 | afatinib | BRD-A58767537 | 0.04 uM",
        "raw_wells": records(trace_raw[["source_file", "original_excel_row",
                                        "Total Cell Count Before Treatment", "Total Cell Count After Treatment",
                                        "Total Control Cell Count", "relative_cell_count", "relative_net_growth"]]),
        "biorepeat_means": records(trace_bio), "equal_weight_mean": records(trace_mean)[0]}
    say(f"  原孔 {len(raw):,} → 生物重复 {len(bio):,} → 库板条件 {len(means):,}；差值 {diffs}")
    say("  示例三批相对计数: " + ", ".join(f"{v:.12f}" for v in trace_bio.relative_cell_count)
        + f"；等权均值 {trace_mean.relative_cell_count.iloc[0]:.12f}")
    # 聚合是观测单位转换而非剔除；漏斗仅从 4140 条库板记录起算。
    funnel(report, 1, "HMS_library_conditions", len(frozen), frozen, "原孔汇总；无库板记录剔除")
    return frozen


def verify_ccle(frozen):
    """克隆内无大 CCLE 源表时，明确记录固定身份/RNA 覆盖表的证据边界。"""
    roots = [REPO_ROOT.parent / "data/raw/CCLE_DepMap_24Q4", DATA / "raw/CCLE_DepMap_24Q4"]
    root = next((p for p in roots if (p / "OmicsExpressionProteinCodingGenesTPMLogp1.csv").is_file()), None)
    info = {"mode": "frozen_reference_coverage_and_baseline", "raw_RNA_verified": False,
            "qualification": "MCF10A 无本发布 RNA 行；不推断其未公开的上游缺失原因"}
    if root is None:
        return info
    path = root / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    header = pd.read_csv(path, nrows=0).columns
    gene_columns = {name.rsplit("(", 1)[-1].rstrip(")"): name for name in header[1:]}
    baseline = pd.read_csv(REF / "baseline.csv", index_col=0)
    selected_columns = [gene_columns[g] for g in baseline.columns]
    pieces = []
    all_ids = []
    for chunk in pd.read_csv(path, usecols=[header[0]] + selected_columns, chunksize=256, index_col=0):
        all_ids.extend(chunk.index)
        pieces.append(chunk.loc[chunk.index.isin(frozen.ModelID.dropna())])
    actual = frozen.ModelID.isin(all_ids)
    assert actual.equals(frozen.actual_RNA_row)
    selected = pd.concat(pieces).reindex(index=baseline.index, columns=selected_columns)
    error = max_diff(selected.to_numpy(), baseline.to_numpy())
    assert error < 1e-12
    model = pd.read_csv(root / "Model.csv")
    missing_ids = frozen.loc[~actual, "ModelID"].unique()
    assert set(missing_ids) <= set(model.ModelID)
    info.update(mode="raw_CCLE_24Q4_RNA", raw_RNA_verified=True, source=str(path),
                RNA_models=len(all_ids), baseline_max_abs_diff=error,
                missing_RNA_models_present_in_Model_csv=missing_ids.tolist())
    return info


def stage2_inputs(frozen, report):
    say("\n[Stage 2/5] 实际 CCLE RNA 行 + 可解析 SMILES")
    ccle = verify_ccle(frozen)
    coverage = pd.read_csv(REF / "drug_structure_coverage.csv")
    pert = raw_table("pert_info").set_index("pert_id")
    assert coverage.pert_id.is_unique
    assert coverage.canonical_smiles.equals(coverage.pert_id.map(pert.canonical_smiles).rename("canonical_smiles"))
    # restricted 不是 SMILES；不把受限结构当作空分子或全零指纹。
    parsed = coverage.canonical_smiles.map(
        lambda value: isinstance(value, str) and value != "restricted" and Chem.MolFromSmiles(value) is not None)
    assert parsed.equals(coverage.parseable_smiles)
    assert frozen.parseable_smiles.equals(frozen.pert_id.map(coverage.set_index("pert_id").parseable_smiles).rename("parseable_smiles"))
    no_rna, no_structure = ~frozen.actual_RNA_row, ~frozen.parseable_smiles
    reason_masks = {"missing_RNA_inclusive": no_rna, "unparseable_SMILES_inclusive": no_structure,
                    "overlap": no_rna & no_structure, "RNA_only": no_rna & ~no_structure,
                    "SMILES_only": ~no_rna & no_structure}
    counts = {name: {"records": int(mask.sum()),
                     "conditions": int(frozen.loc[mask, "biological_condition_group"].nunique())}
              for name, mask in reason_masks.items()}
    kept = frozen[frozen.actual_RNA_row & frozen.parseable_smiles].copy()
    assert len(kept) == 3420
    dropped_drugs = records(coverage[~parsed])
    report["stage2"] = {"drop_counts": counts, "dropped_drugs": dropped_drugs,
                        "missing_RNA_cells": records(frozen.loc[no_rna, ["cell_id", "ModelID"]].drop_duplicates()),
                        "CCLE_verification": ccle, "parseability_recomputed_with_rdkit": True,
                        "sequential_drops": {"missing_RNA": int(no_rna.sum()),
                                             "SMILES_after_RNA_filter": int((no_structure & ~no_rna).sum())}}
    say(f"  RNA缺行（含重叠）{counts['missing_RNA_inclusive']}，结构受限（含重叠）"
        f"{counts['unparseable_SMILES_inclusive']}，重叠 {counts['overlap']}")
    say(f"  剔除药物: {dropped_drugs}")
    funnel(report, 2, "CCLE_background_and_parseable_structure", len(frozen), kept,
           "MCF10A 缺 RNA 690条；其余细胞 WH-4-025/restricted 30条；共720条")
    return kept


def extract_gctx(path, names, genes, label):
    """忠实采用冻结脚本的 HDF5 方向及 256 行提取，再按 landmark 顺序取列。"""
    with h5py.File(path, "r") as f:
        row_ids = f["0/META/COL/id"].asstr()[:]
        gene_ids = f["0/META/ROW/id"].asstr()[:]
        row_map = {value: i for i, value in enumerate(row_ids)}
        gene_map = {value: i for i, value in enumerate(gene_ids)}
        assert len(row_map) == len(row_ids) and len(gene_map) == len(gene_ids)
        cols = np.array([gene_map[g] for g in genes])
        order = sorted(range(len(names)), key=lambda i: row_map[names[i]])
        assert len(names) == len(set(names))
        result = np.empty((len(names), len(genes)), dtype=np.float32)
        matrix = f["0/DATA/0/matrix"]
        for start in range(0, len(order), 256):
            selection = order[start:start + 256]
            block = matrix[[row_map[names[i]] for i in selection], :]
            result[selection] = block[:, cols]
            if start % 2048 == 0 or start + 256 >= len(order):
                say(f"  {label}: 已读取 {min(start + 256, len(order)):,}/{len(order):,} 行 × 978基因")
        metadata = {"source": str(path), "matrix_shape": list(matrix.shape),
                    "rows_extracted": len(names), "genes_extracted": len(genes),
                    "HDF5_direction": "META/COL=id 为 matrix 行；META/ROW=id 为基因列",
                    "file_bytes": path.stat().st_size}
    return result, metadata


def stage3_expression(c, unpacked, skip, report):
    say("\n[Stage 3/5] 同库板 24 h Level 5 签名与 978 landmark 基因对齐")
    before = len(c)
    missing = c[c.signature_ids_24h.isna()]
    # 追查缺行到原始签名表，区分本发布无对应签名与提取脚本漏读。
    raw_sig = raw_table("sig_info")
    raw_sig["library_id"] = raw_sig.sig_id.str.split("_").str[0]
    raw_sig["dose_uM"] = pd.to_numeric(raw_sig.pert_idose.str.split().str[0], errors="coerce")
    raw_sig["time_h"] = pd.to_numeric(raw_sig.pert_itime.str.split().str[0], errors="coerce")
    absent = []
    for row in missing.itertuples():
        same = raw_sig[(raw_sig.library_id.eq(row.library_id)) & raw_sig.cell_id.eq(row.cell_id)
                       & raw_sig.pert_id.eq(row.pert_id) & np.isclose(raw_sig.dose_uM, row.dose_uM)]
        assert not same.time_h.eq(24).any()
        item = {key: getattr(row, key) for key in KEYS + ["pert_iname", "library_condition_id"]}
        item.update(available_times_h=same.time_h.tolist(), available_signature_ids=same.sig_id.tolist(),
                    reason="原始 sig_info 在该库板/细胞/药物/剂量下无24h签名；上游未公开原因未推断")
        absent.append(item)
    c = c[c.signature_ids_24h.notna()].copy().rename(columns={"signature_ids_24h": "sig_id"})
    signatures = pd.read_csv(REF / "lincs_signatures.csv")
    c = c.merge(signatures, on=["sig_id", "cell_id", "pert_id", "dose_uM"], validate="one_to_one")
    assert len(c) == len(signatures) == 3410 and c.time_h.eq(24).all()
    assert c.sig_id.str.split("_").str[0].eq(c.library_id).all()
    source_sig = raw_sig.set_index("sig_id").loc[c.sig_id]
    for key in ["cell_id", "pert_id", "distil_id"]:
        assert np.array_equal(source_sig[key], c[key])
    assert np.array_equal(source_sig.dose_uM, c.dose_uM) and source_sig.time_h.eq(24).all()
    expression = npz(REF / "lincs_expression.npz")
    features = npz(REF / "input_features.npz")
    annotation = raw_table("gene_info")
    genes = annotation.loc[annotation.pr_is_lm.eq(1), "pr_gene_id"].to_numpy(dtype=str)
    assert len(genes) == 978
    assert np.array_equal(genes, expression["gene_id"]) and np.array_equal(genes, features["gene_id"])
    assert np.array_equal(expression["sig_id"], signatures.sig_id.to_numpy(dtype=str))
    assert np.isfinite(expression["expression_z"]).all()
    paths = sorted((unpacked / "LINCS_GSE70138").glob("*Level5*.gctx"))
    l5 = {"mode": "archived_reference_expression", "raw_verified": False,
          "reason": "--skip-l3-gctx 同时跳过大 L3/L5 矩阵" if skip else "未发现解压的 Level5 GCTX"}
    if paths and not skip:
        path = only_file(unpacked / "LINCS_GSE70138", "*Level5*.gctx")
        rebuilt, metadata = extract_gctx(path, signatures.sig_id.tolist(), genes, "Level5")
        assert np.array_equal(rebuilt, expression["expression_z"])
        l5 = dict(metadata, mode="full_raw_l5_gctx", raw_verified=True, max_abs_diff=0.0)
        expression["expression_z"] = rebuilt
    report["stage3"] = {"missing_24h_candidates": absent, "landmark_order_exact": True,
                        "landmark_genes": len(genes), "signature_metadata_exact": True, "L5": l5}
    funnel(report, 3, "24h_signature_available", before, c, "原始sig_info无匹配24h签名：10条（明细见JSON）")
    return c, expression, features


def stage4_qc(c, expression, unpacked, skip, report):
    say("\n[Stage 4/5] Level 3 逐实例有效性与同检测板 DMSO 对照 QC")
    inst = raw_table("inst_info")
    assert inst.inst_id.is_unique
    refs = c[["sig_id", "cell_id", "pert_id", "dose_uM", "distil_id"]].assign(
        inst_id=lambda d: d.distil_id.str.split("|")).explode("inst_id")
    before_join = len(refs)
    refs = refs.merge(inst, on="inst_id", suffixes=("", "_instance"), validate="many_to_one")
    assert len(refs) == before_join == 10512 and refs.inst_id.nunique() == 10512
    assert refs.cell_id.eq(refs.cell_id_instance).all() and refs.pert_id.eq(refs.pert_id_instance).all()
    assert np.isclose(refs.dose_uM, refs.pert_dose).all() and refs.pert_time.eq(24).all()
    assert refs.pert_dose_unit.str.lower().eq("um").all() and refs.pert_type.eq("trt_cp").all()
    assert refs.pert_time_unit.eq("h").all()
    controls = inst[inst.pert_type.eq("ctl_vehicle") & inst.pert_time.eq(24)
                    & inst.det_plate.isin(refs.det_plate)].copy()
    assert len(controls) == 799 and controls.det_plate.nunique() == refs.det_plate.nunique() == 32
    assert controls.pert_id.eq("DMSO").all() and controls.pert_iname.eq("DMSO").all()
    archived_refs = pd.read_csv(TRAIN / "instance_checks.csv")
    archived_ctrl = pd.read_csv(TRAIN / "vehicle_plate_checks.csv")
    archived_c = pd.read_csv(TRAIN / "all_24h_candidates_with_QC.csv")
    # 先完整核对实例来源及顺序，再允许轻量模式使用归档的数值有效标志。
    exact_frame(refs, archived_refs[refs.columns], "实例原始元数据")
    ctrl_counts = controls.groupby("det_plate").size()
    assert np.array_equal(ctrl_counts.index, archived_ctrl.det_plate)
    assert np.array_equal(ctrl_counts, archived_ctrl.listed_vehicle_instances)
    paths = sorted((unpacked / "LINCS_GSE70138").glob("*Level3*.gctx"))
    full = bool(paths) and not skip
    if not full and not skip:
        raise FileNotFoundError(f"未找到 {unpacked}/LINCS_GSE70138/*Level3*.gctx；"
                                "请提供 --unpacked-root 或显式使用 --skip-l3-gctx。")
    report["l3_verification_mode"] = "full_raw_l3_gctx" if full else "archived_qc_tables"
    extracted = read_map = medians = None
    if full:
        needed = sorted(set(refs.inst_id) | set(controls.inst_id))
        assert len(needed) == 11311
        extracted, l3meta = extract_gctx(only_file(unpacked / "LINCS_GSE70138", "*Level3*.gctx"),
                                       needed, expression["gene_id"], "Level3")
        l3meta["raw_verified"] = True
        read_map = {name: i for i, name in enumerate(needed)}
        ctrl_info, medians = [], {}
        for plate, group in controls.groupby("det_plate"):
            a = extracted[[read_map[name] for name in group.inst_id]]
            valid = np.isfinite(a).all(1) & (a != -666).all(1)
            good = a[valid]
            if len(good):
                medians[plate] = np.median(good, axis=0)
            ctrl_info.append({"det_plate": plate, "listed_vehicle_instances": len(a),
                              "finite_non_sentinel_vehicle_instances": int(valid.sum()),
                              "vehicle_median_gene_sd": float(np.median(np.std(good, axis=0))) if len(good) > 1 else np.nan})
        ctrl = pd.DataFrame(ctrl_info)
        refs["L3_finite_non_sentinel"] = [bool(np.isfinite(extracted[read_map[name]]).all()
                                                 and (extracted[read_map[name]] != -666).all()) for name in refs.inst_id]
        l3meta["vehicle_sd_max_abs_diff_vs_archived"] = max_diff(ctrl.vehicle_median_gene_sd,
                                                                archived_ctrl.vehicle_median_gene_sd)
        assert l3meta["vehicle_sd_max_abs_diff_vs_archived"] < 1e-12
        for key in ["det_plate", "listed_vehicle_instances", "finite_non_sentinel_vehicle_instances"]:
            assert np.array_equal(ctrl[key], archived_ctrl[key])
    else:
        ctrl = archived_ctrl.copy()
        refs["L3_finite_non_sentinel"] = archived_refs.L3_finite_non_sentinel.to_numpy()
        l3meta = {"raw_verified": False,
                  "qualification": "从3张归档QC表交叉核对；未重新测量原始L3数值或重算相关系数"}
    ctrl_by_plate = ctrl.set_index("det_plate")
    assert ctrl.finite_non_sentinel_vehicle_instances.ge(2).all()
    assert ctrl.finite_non_sentinel_vehicle_instances.le(ctrl.listed_vehicle_instances).all()
    refs["valid_vehicle_count_same_plate"] = refs.det_plate.map(ctrl_by_plate.finite_non_sentinel_vehicle_instances)
    refs["instance_technical_valid"] = refs.L3_finite_non_sentinel & refs.valid_vehicle_count_same_plate.ge(2)
    exact_frame(refs, archived_refs, "逐实例QC")
    quality = []
    for sid, group in refs.groupby("sig_id"):
        good = group[group.instance_technical_valid]
        item = {"sig_id": sid, "n_listed_instances": len(group),
                "n_technical_valid_instances": len(good),
                "all_instances_technical_valid": bool(group.instance_technical_valid.all())}
        if full:
            delta = np.array([extracted[read_map[r.inst_id]] - medians[r.det_plate] for r in good.itertuples()])
            correlations = []
            for i in range(len(delta)):
                for j in range(i + 1, len(delta)):
                    if np.std(delta[i]) > 0 and np.std(delta[j]) > 0:
                        correlations.append(float(np.corrcoef(delta[i], delta[j])[0, 1]))
            item["replicate_vehicle_centered_pearson_median"] = float(np.median(correlations)) if correlations else np.nan
        quality.append(item)
    q = pd.DataFrame(quality).set_index("sig_id").reindex(c.sig_id)
    archived_q = archived_c.set_index("sig_id").reindex(c.sig_id)
    for key in ["n_listed_instances", "n_technical_valid_instances", "all_instances_technical_valid"]:
        assert np.array_equal(q[key], archived_q[key])
    if full:
        error = max_diff(q.replicate_vehicle_centered_pearson_median,
                         archived_q.replicate_vehicle_centered_pearson_median)
        assert error < 1e-12
        l3meta["pearson_median_max_abs_diff_vs_archived"] = error
    # 相关系数只作诊断，未用于筛选。核对重算值后，保留CSV冻结的小数精度，
    # 避免读写十进制带来的 1e-16 差异破坏 samples.csv 的逐值精确对照。
    q["replicate_vehicle_centered_pearson_median"] = archived_q.replicate_vehicle_centered_pearson_median
    c = c.merge(q.reset_index()[["sig_id"] + QC_COLUMNS], on="sig_id", validate="one_to_one")
    assert np.array_equal(c.distil_nsample, c.n_listed_instances)
    c["main_eligible"] = c.all_instances_technical_valid & c.n_technical_valid_instances.ge(2)
    c["eligibility_reason"] = np.where(c.main_eligible,
        "source_and_numeric_checks_passed;at_least_two_instances_with_same_plate_controls",
        "fewer_than_two_valid_instances_or_control_issue")
    candidate_parity = exact_frame(c, archived_c, "全体3410候选及QC")
    excluded = c[~c.main_eligible]
    assert len(excluded) == 46 and excluded.n_listed_instances.eq(1).all()
    assert c.all_instances_technical_valid.all() and refs.L3_finite_non_sentinel.all()
    assert set(excluded.sig_id) == set(archived_c.loc[~archived_c.main_eligible, "sig_id"])
    kept = c[c.main_eligible].sort_values("library_condition_id").reset_index(drop=True)
    assert (len(kept), kept.biological_condition_group.nunique(), kept.cell_id.nunique(), kept.pert_id.nunique()) == (3364, 3133, 5, 106)
    report["stage4"] = {"L3": l3meta, "treatment_instances": len(refs), "DMSO_controls": len(controls),
                        "det_plates": len(ctrl), "all_treatments_finite_non_sentinel": True,
                        "numeric_scope": "每个实例978个实测landmark；不对11350个推断基因做QC",
                        "eligibility_rule": "全部列出实例技术有效 AND 至少2个有效处理实例；每实例同板有效vehicle至少2个",
                        "instance_identity_dose_time_unit_checks": True,
                        "all_vehicles_finite_non_sentinel": bool(ctrl.listed_vehicle_instances.eq(ctrl.finite_non_sentinel_vehicle_instances).all()),
                        "min_valid_same_plate_controls": int(ctrl.finite_non_sentinel_vehicle_instances.min()),
                        "max_valid_same_plate_controls": int(ctrl.finite_non_sentinel_vehicle_instances.max()),
                        "eligible_replicate_count_distribution": kept.n_technical_valid_instances.value_counts().sort_index().to_dict(),
                        "excluded_single_instance_signatures": records(excluded[["library_condition_id", "biological_condition_group",
                            "sig_id", "library_id", "cell_id", "pert_id", "pert_iname", "dose_uM",
                            "n_listed_instances", "n_technical_valid_instances", "eligibility_reason"]]),
                        "candidate_table_parity": candidate_parity,
                        "diagnostic_precision": "相关系数在全量模式重算核对<1e-12后采用归档CSV精度；不参与入组",
                        "no_TAS_or_expression_magnitude_filter": True}
    report["example_trace"]["L3_instances"] = records(refs[refs.sig_id.eq(kept.loc[example_mask(kept), "sig_id"].iloc[0])])
    funnel(report, 4, "technical_QC_at_least_two_instances", len(c), kept, "46条单实例签名；无数值无效或同板对照不足导致的额外剔除")
    say(f"  模式={report['l3_verification_mode']}；处理实例10,512；DMSO799；检测板32；单实例剔除46")
    return kept


def assign_folds(frame, key):
    """与 src/prepare.py 相同：按不同生物条件负载分组，固定种子哈希破同分。"""
    counts = frame.drop_duplicates("biological_condition_group").groupby(key).size()
    loads, assignment = [0] * 5, {}
    for group, size in sorted(counts.items(), key=lambda item: (
            -item[1], hashlib.sha256(f"20261107|{item[0]}".encode()).hexdigest())):
        fold = min(range(5), key=lambda j: (loads[j], j))
        assignment[group] = fold
        loads[fold] += int(size)
    return frame[key].map(assignment)


def make_splits(c):
    parts = []
    for protocol, key in [("unseen_drug", "chemical_split_group"), ("unseen_scaffold", "scaffold_split_group")]:
        assignment = assign_folds(c, key)
        for fold in range(5):
            d = c[["sample_index", "library_condition_id", "biological_condition_group", key]].rename(columns={key: "heldout_group"}).copy()
            d["protocol"], d["fold"] = protocol, fold
            d["role"] = np.where(assignment.eq(fold), "test",
                                 np.where(assignment.eq((fold + 1) % 5), "validation", "train"))
            parts.append(d)
    d = c[["sample_index", "library_condition_id", "biological_condition_group"]].copy()
    d["heldout_group"] = c.biological_condition_group
    d["protocol"], d["fold"] = "dose_extrapolation", 0
    assert c.dose_uM.isin([0.04, 0.12, 0.37, 1.11, 3.33, 10.0]).all()
    d["role"] = np.where(np.isclose(c.dose_uM, 10), "test",
                         np.where(np.isclose(c.dose_uM, 3.33), "validation", "train"))
    assert c.loc[d.role.eq("train"), "dose_uM"].le(1.11).all()
    parts.append(d)
    return pd.concat(parts, ignore_index=True)


def check_leakage_and_scaling(c, splits, arrays):
    checks, scaling = [], []
    for (protocol, fold), d in splits.groupby(["protocol", "fold"]):
        assert len(d) == len(c) and d.sample_index.is_unique
        assert set(d.sample_index) == set(c.sample_index)
        assert set(d.role) == {"train", "validation", "test"}
        keys = ["biological_condition_group"]
        if protocol in ("unseen_drug", "unseen_scaffold"):
            keys.append("heldout_group")
        for key in keys:
            leaking = int(d.groupby(key).role.nunique().gt(1).sum())
            assert leaking == 0
            checks.append({"protocol": protocol, "fold": int(fold), "group": key, "leaking_groups": leaking})
        # 标准化不写回训练数组；按 train.py 的规则逐折演示并检查参数。
        train = d.loc[d.role.eq("train"), "sample_index"].to_numpy()
        w = arrays["condition_weight"][train]
        cell = arrays["baseline"][arrays["cell_index"]].astype(np.float32)
        cm = np.average(cell[train], axis=0, weights=w)
        cs = np.sqrt(np.average((cell[train] - cm) ** 2, axis=0, weights=w)).clip(1e-4)
        dm = float(np.average(arrays["log10_dose_uM"][train], weights=w))
        ds = float(np.sqrt(np.average((arrays["log10_dose_uM"][train] - dm) ** 2, weights=w)))
        fm = float(np.average(arrays["count72"][train], weights=w))
        fs = float(np.sqrt(np.average((arrays["count72"][train] - fm) ** 2, weights=w)))
        assert np.isfinite(cm).all() and np.isfinite(cs).all() and ds > 0 and fs > 0
        scaling.append({"protocol": protocol, "fold": int(fold), "fit_role": "train",
                        "fit_records": len(train), "fit_weight_sum": float(w.sum()),
                        "cell_mean_range": [float(cm.min()), float(cm.max())],
                        "cell_std_range": [float(cs.min()), float(cs.max())],
                        "cell_mean_sha256": hashlib.sha256(cm.tobytes()).hexdigest(),
                        "cell_std_sha256": hashlib.sha256(cs.tobytes()).hexdigest(),
                        "dose_mean": dm, "dose_std": ds, "count72_mean": fm, "count72_std": fs})
    return checks, scaling


def stage5_assembly(c, expression, features, write_dir, report):
    say("\n[Stage 5/5] 1/n_plates 权重、张量、确定性划分、逐值 parity 与零泄漏")
    drugs, cells = pd.read_csv(REF / "drugs.csv"), pd.read_csv(REF / "cells.csv")
    c = c.merge(drugs[["pert_id", "drug_index", "chemical_split_group", "scaffold_split_group"]],
                on="pert_id", validate="many_to_one")
    c = c.merge(cells[["cell_id", "ModelID", "cell_index"]], on=["cell_id", "ModelID"], validate="many_to_one")
    c["sample_index"] = np.arange(len(c))
    c["condition_weight"] = 1.0 / c.groupby("biological_condition_group").library_condition_id.transform("size")
    c["sensitivity_rep3plus"] = c.n_technical_valid_instances.ge(3)
    c["sensitivity_without_CGP60474"] = c.pert_id.ne("BRD-K79090631")
    assert np.array_equal(cells.cell_index, np.arange(5)) and np.array_equal(drugs.drug_index, np.arange(106))
    assert c.groupby("biological_condition_group").condition_weight.sum().eq(1.0).all()
    assert c.condition_weight.sum() == 3133.0
    sizes = c.groupby("biological_condition_group").size()
    assert sizes.isin([1, 2]).all()
    # 指纹采用已冻结的标准化母体结构；不更改 v0 的化学/骨架分组定义。
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    morgan = np.stack([generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))
                       for s in drugs.standardized_isomeric_smiles])
    assert np.array_equal(morgan, features["drug_morgan"]) and set(np.unique(morgan)) == {0, 1}
    zmap = {sid: i for i, sid in enumerate(expression["sig_id"])}
    genes = expression["gene_id"].astype(str)
    baseline = pd.read_csv(REF / "baseline.csv", index_col=0).reindex(index=cells.ModelID, columns=genes)
    arrays = dict(expression_z24=expression["expression_z"][[zmap[sid] for sid in c.sig_id]],
                  baseline=baseline.to_numpy(dtype=np.float32), drug_morgan=morgan,
                  cell_index=c.cell_index.to_numpy(), drug_index=c.drug_index.to_numpy(),
                  log10_dose_uM=np.log10(c.dose_uM.to_numpy()).astype(np.float32),
                  count72=c.relative_cell_count.to_numpy(), condition_weight=c.condition_weight.to_numpy(),
                  gene_id=genes, library_condition_id=c.library_condition_id.to_numpy(dtype=str))
    frozen = npz(TRAIN / "arrays.npz")
    assert set(arrays) == set(frozen)
    parity, stats = {}, {}
    units = {"baseline": "log2(TPM+1)", "drug_morgan": "binary bit",
             "log10_dose_uM": "log10(dose / 1 uM)", "expression_z24": "L1000 Level5 z-score, 24 h",
             "count72": "Nt/Nc, 72 h, dimensionless", "condition_weight": "1/n_library_plates",
             "cell_index": "zero-based index", "drug_index": "zero-based index",
             "gene_id": "Entrez gene ID", "library_condition_id": "HMS library condition ID"}
    for key, array in arrays.items():
        assert array.dtype == frozen[key].dtype and np.array_equal(array, frozen[key]), key
        parity[key] = {"exact": True, "max_abs_diff": 0.0 if array.dtype.kind in "fiu" else None}
        stats[key] = {"shape": list(array.shape), "dtype": str(array.dtype), "unit": units[key]}
        if array.dtype.kind in "fiu":
            assert np.isfinite(array).all()
            stats[key].update(min=float(array.min()), max=float(array.max()),
                              mean=float(array.mean()), std=float(array.std()), finite=True)
    assert arrays["baseline"].shape == (5, 978) and arrays["drug_morgan"].shape == (106, 2048)
    assert arrays["expression_z24"].shape == (3364, 978)
    samples_parity = exact_frame(c, pd.read_csv(TRAIN / "samples.csv"), "samples.csv")
    splits = make_splits(c)
    splits_parity = exact_frame(splits, pd.read_csv(TRAIN / "splits.csv"), "splits.csv")
    leakage, scaling = check_leakage_and_scaling(c, splits, arrays)
    joined = splits.merge(c[["sample_index", "pert_id", "cell_id"]], on="sample_index", validate="many_to_one")
    summary = joined.groupby(["protocol", "fold", "role"]).agg(
        records=("sample_index", "size"), conditions=("biological_condition_group", "nunique"),
        n_drugs=("pert_id", "nunique"), n_cells=("cell_id", "nunique")).reset_index()
    report["stage5"] = {"weights": {"sum": float(c.condition_weight.sum()),
                         "one_plate_conditions": int(sizes.eq(1).sum()), "one_plate_records": int(c.condition_weight.eq(1).sum()),
                         "two_plate_conditions": int(sizes.eq(2).sum()), "two_plate_records": int(c.condition_weight.eq(0.5).sum())},
                        "feature_stats": stats, "arrays_parity": parity, "samples_parity": samples_parity,
                        "splits_parity": splits_parity, "leakage_checks": leakage,
                        "split_rule": "不同生物条件数降序；sha256('20261107|group')升序破同分；分给(load,fold)最小折；test=f, validation=(f+1)%5",
                        "dose_split_rule_uM": {"train": "<=1.11", "validation": 3.33, "test": 10.0},
                        "scaling_rule": {"fit_role": "train only", "weights": "condition_weight",
                            "formula": "mu=sum(w*x)/sum(w); sigma=sqrt(sum(w*(x-mu)^2)/sum(w)); (x-mu)/sigma",
                            "baseline": "按训练样本展开到978维；每基因加权标准化，std最小1e-4",
                            "dose": "log10_dose_uM训练折加权均值/标准差",
                            "count72": "训练折加权均值/标准差；预测后逆变换",
                            "Morgan": "保持0/1，不标准化", "expression_z24": "保持已发布Level5 z-score，不二次标准化",
                            "arrays_on_disk": "保存原始尺度；标准化在每折训练时进行"},
                        "train_scaling_by_fold": scaling, "Morgan_recomputed_exact": True}
    one = c[example_mask(c)]
    assert len(one) == 1
    sample = one.iloc[0]
    index = int(sample.sample_index)
    report["example_trace"].update(final_sample=records(one)[0],
        features={"sample_index": index, "cell_index": int(sample.cell_index), "drug_index": int(sample.drug_index),
                  "log10_dose_uM": float(arrays["log10_dose_uM"][index]),
                  "Morgan_on_bits": int(arrays["drug_morgan"][sample.drug_index].sum()),
                  "baseline_first5": arrays["baseline"][sample.cell_index, :5].tolist(),
                  "gene_id_first5": genes[:5].tolist(),
                  "expression_z24_first5": arrays["expression_z24"][index, :5].tolist(),
                  "count72": float(arrays["count72"][index]), "condition_weight": float(sample.condition_weight)},
        split_roles=records(splits[splits.sample_index.eq(index)][["protocol", "fold", "role", "heldout_group"]]))
    funnel(report, 5, "weights_features_splits_verified", len(c), c, "无剔除；精确parity与零组泄漏核对")
    if write_dir is not None:
        write_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(write_dir / "arrays.npz", **arrays)
        c.to_csv(write_dir / "samples.csv", index=False)
        splits.to_csv(write_dir / "splits.csv", index=False)
        saved = npz(write_dir / "arrays.npz")
        assert all(saved[key].dtype == array.dtype and np.array_equal(saved[key], array)
                   for key, array in arrays.items())
        exact_frame(pd.read_csv(write_dir / "samples.csv"), c, "另存samples的读回核对")
        exact_frame(pd.read_csv(write_dir / "splits.csv"), splits, "另存splits的读回核对")
        report["rebuilt_dir"] = str(write_dir)
        report["rebuilt_files_readback_exact"] = True
    else:
        report["rebuilt_dir"] = None
    say(f"  数组10项、samples全部{len(c.columns)}列、splits全部7列：逐值100%一致")
    say(f"  权重和3133.0；单板{sizes.eq(1).sum()}条件；双板{sizes.eq(2).sum()}条件；11套划分零组泄漏")
    say("  标准化仅拟合train：baseline逐基因/log10剂量/count72按condition_weight加权；Morgan保持0/1")
    return summary


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unpacked-root", type=Path, help="包含LINCS_GSE70138子目录的解压根目录")
    parser.add_argument("--skip-l3-gctx", action="store_true", help="轻量模式：跳过L3/L5大矩阵，交叉核对归档QC表")
    parser.add_argument("--write-rebuilt-dir", type=Path, help="可选：另存arrays.npz/samples.csv/splits.csv；不覆盖冻结输入")
    args = parser.parse_args()
    roots = [REPO_ROOT.parent / "data/unpacked", DATA / "unpacked"]
    args.unpacked_root = (args.unpacked_root or next((p for p in roots if p.is_dir()), roots[0])).resolve()
    if args.write_rebuilt_dir is not None:
        args.write_rebuilt_dir = args.write_rebuilt_dir.resolve()
        # 输出不能覆盖任何冻结目录；仓库内仅允许写本迭代。
        if (args.write_rebuilt_dir.is_relative_to(REPO_ROOT)
                and not args.write_rebuilt_dir.is_relative_to(HERE)):
            parser.error("仓库内 --write-rebuilt-dir 只能指向 iterations/01_v0_preprocessing/ 下的新目录")
        if args.write_rebuilt_dir.is_relative_to(REF) or args.write_rebuilt_dir.is_relative_to(TRAIN):
            parser.error("不可覆盖 data/reference 或 data/training")
    return args


def main():
    args = parse_args()
    report = {"iteration": "iter-01", "pipeline": "initial_v0_frozen_preprocessing",
              "repo_branch_expected": "step-by-step-repro", "unpacked_root": str(args.unpacked_root),
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "pandas": pd.__version__, "h5py": h5py.__version__,
                              "rdkit": rdBase.rdkitVersion, "xlrd": importlib.metadata.version("xlrd") if importlib.util.find_spec("xlrd") else None},
              "frozen_definition_boundary": "身份、化学/骨架分组、候选签名映射采用已审计data/reference；不扩展或重定义v0",
              "funnel": []}
    say(f"v0数据预处理审计；解压目录={args.unpacked_root}；轻量模式={args.skip_l3_gctx}")
    c = stage1_hms(report)
    c = stage2_inputs(c, report)
    c, expression, features = stage3_expression(c, args.unpacked_root, args.skip_l3_gctx, report)
    c = stage4_qc(c, expression, args.unpacked_root, args.skip_l3_gctx, report)
    summary = stage5_assembly(c, expression, features, args.write_rebuilt_dir, report)
    # 仅在全部断言通过后发布本次审计表；不运行会改写data/的冻结脚本。
    inputs = sorted(REF.glob("*")) + [TRAIN / name for name in ["arrays.npz", "samples.csv", "splits.csv",
              "all_24h_candidates_with_QC.csv", "instance_checks.csv", "vehicle_plate_checks.csv"]]
    inputs += sorted((DATA / "raw/LINCS_GSE70138").glob("*.gz"))
    inputs += [REPO_ROOT / "src" / name for name in ["aggregate_hms.py", "extract_reference.py", "prepare.py", "check.py", "train.py"]]
    report["input_sha256"] = {str(p.relative_to(REPO_ROOT)): digest(p) for p in inputs if p.is_file()}
    report["success"] = True
    report["final_cohort"] = {"records": len(c), "conditions": int(c.biological_condition_group.nunique()),
                              "cells": sorted(c.cell_id.unique()), "n_cells": int(c.cell_id.nunique()),
                              "n_drugs": int(c.pert_id.nunique()), "expression_time_h": 24, "count_time_h": 72}
    pd.DataFrame(report["funnel"]).to_csv(HERE / "v0_cohort_funnel.csv", index=False)
    summary.to_csv(HERE / "v0_split_summary.csv", index=False)
    (HERE / "v0_preprocessing_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    say("\n全部核对通过。已写入本迭代的3份小型审计文件；未训练模型。")


if __name__ == "__main__":
    main()
