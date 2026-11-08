#!/usr/bin/env python3
"""Create an overview and unified mapping table from audited source outputs."""
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "audit"
d = json.loads((A / "audit_results.json").read_text())
m = json.loads((A / "mixseq_pilot_audit.json").read_text())
h = d["HMS_LINCS"]
i = d["HMS_input_metadata_availability"]
records = []
for row in pd.read_csv(A / "hms_cell_mapping.csv").fillna("").to_dict("records"):
    records.append(dict(entity_type="cell", source="HMS", source_id=row["cell_id"],
                        canonical_id=row["RRID"], canonical_id_type="Cellosaurus",
                        depmap_id=row["ModelID"], label=row["CellLineName"],
                        structure_smiles="", rna_profile_id=row["ProfileID"],
                        mapping_basis="Pinned_DepMap_24Q4_Model_and_default_profile_tables",
                        unresolved="" if row["ProfileID"] else "No_default_RNA_profile_in_this_pinned_release"))
for row in pd.read_csv(A / "hms_drug_mapping.csv").fillna("").to_dict("records"):
    smi = str(row.get("canonical_smiles", ""))
    if smi in ["-666", "nan"]:
        smi = ""
    records.append(dict(entity_type="compound", source="HMS", source_id=row["pert_id"],
                        canonical_id=row["pert_id"], canonical_id_type="Broad_compound_ID",
                        depmap_id="", label=row["drug_name"], structure_smiles=smi,
                        rna_profile_id="", mapping_basis="HMS_Broad_ID_join_to_GSE92742_official_pert_info",
                        unresolved="SMILES_not_chemically_validated" if smi else "Structure_missing_in_downloaded_metadata"))
pd.DataFrame(records).to_csv(ROOT / "data_mapping.csv", index=False, encoding="utf-8-sig")
report = f"""# 实际匹配规模与可执行任务

统计日期：2026年11月5日。所有数字由本地脚本重算。候选元数据对应、矩阵可读和可训练样本是不同阶段。

## Broad HMS 与 LINCS

| 项目 | 数量 | 解释 |
|---|---:|---|
| 功能原始表行数 | {h['raw_observations']:,} | 含技术重复与重复条件 |
| 细胞系 | {h['cells']} | 含非肿瘤MCF10A |
| 药物Broad ID | {h['compounds']} | 不代表全部已有结构输入 |
| 独立细胞药物组合 | {h['matched_cell_drug_pairs']} | 与剂量条件数区分 |
| 唯一功能剂量条件 | {h['unique_functional_conditions']:,} | 全部在LJP表达元数据中找到条件对应 |
| 表达签名对应记录 | {h['matched_LJP_signature_records']:,} | 跨来源与时间，不能当独立功能样本 |
| 同版本3小时和24小时条件 | {h['matched_3_to_24_conditions']:,} | GSE70138 |
| 同版本6小时和24小时条件 | {h['matched_6_to_24_conditions']:,} | GSE92742；不与上行直接相加 |
| 当前已取得结构字符串的药物 | {i['compounds_with_structure_string_in_downloaded_metadata']} | 缺{i['compounds_without_structure_string']}个；字符串还未做化学有效性检查 |
| 有24Q4默认RNA profile索引的细胞 | {i['cells_with_default_rna_profile_index']} | 六细胞均有模型ID；MCF10A无该索引 |
| 同时具备上述结构与RNA索引的功能条件 | {i['functional_conditions_with_both_input_metadata']:,} | 还不是矩阵QC后的可训练数 |
| 同时具备上述输入及3/24小时记录的条件 | {i['within_release_3_to_24_conditions_with_both_input_metadata']:,} | 结构补齐或更换基线来源后可能改变 |

表达端矩阵、控制与批次来源未完成核验，因此不报告完整联合预测训练样本数。

## MIX seq 原始表达试读

| 项目 | 数量或内容 |
|---|---|
| 药物 | Trametinib |
| 表达剂量 | {m['trametinib_expression_dose_uM']} μM；来自原始补充表 |
| 时间 | 6小时和24小时，分别有DMSO对照 |
| normal单细胞合计 | {m['qc_normal_total']:,} |
| 四个条件各至少20个normal细胞的细胞系 | {m['all_four_conditions_ge_20_cells']} |
| 上述细胞系且有作者敏感性标签 | {m['ge_20_and_sensitivity_label']} |
| 功能标签 | 1−归一化AUC，属于曲线敏感性 |

可以开展条件级早期表达、异质性与敏感性关联验证。不能把这些细胞当作9,894个独立功能样本，也不能把AUC解释成0.1 μM指定终点的活性。

## 其余匹配任务

| 任务 | 当前状态 |
|---|---|
| LINCS–CTRP浓度级匹配 | 文献和作者代码已核对；原始响应表本轮未恢复；不报新交集数 |
| LINCS–GDSC | 官方raw入口已定位；尚未读取全部响应表 |
| LINCS–PRISM | 两张元数据和转换说明已读；响应矩阵及非缺失交集未统计 |
| Tahoe–药敏资源 | 已读三张注释表；细胞注释范围大于实际表达面板；不能据此计算实际交集 |
| sci-Plex独立功能标签 | 原文回收数代理与独立验证已区分；完整验证覆盖仍待读补充表 |
| NCI时间表达–功能 | 原始入口和设计已核验；文件级匹配未进行 |

## 对应文件

- [HMS候选配对](audit/hms_lincs_candidate_pairs.csv)
- [具备当前输入元数据的HMS候选配对](audit/hms_complete_input_metadata_candidate_pairs.csv)
- [MIX-seq表达与敏感性对应](audit/mixseq_pilot_expression_sensitivity_pairs.csv)
- [统一细胞与药物映射](data_mapping.csv)
- [原始下载地址](download_sources.csv)
- [文件校验清单](audit/downloaded_file_manifest.csv)

原始来源：[HMS论文](https://doi.org/10.1038/s41467-017-01383-w)、[HMS功能表目录](https://lincs.hms.harvard.edu/db/datasets/)、[LINCS Phase I](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)、[LINCS Phase II](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)、[MIX-seq原始数据](https://doi.org/10.6084/m9.figshare.10298696)、[DepMap 24Q4文件索引](https://api.figshare.com/v2/articles/27993248)。
"""
(ROOT / "overlap_report.md").write_text(report, encoding="utf-8")
print("Saved overlap_report.md and", len(records), "mapping rows")
