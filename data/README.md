# 数据文件：原始观测、对接依据与训练集

先读[每个数据集是什么、为什么使用](../docs/02_data.md)，再读[实际对接过程](../docs/03_integration.md)。

## training：实际送入模型的数据

| 文件/数组字段 | 数量 | 含义 |
|---|---|---|
| `samples.csv` | 3,364行 | 库板观测的条件、功能答案、表达签名、QC与来源；`sample_index`对应数组行号 |
| `splits.csv` | 按协议与折展开 | `role`区分训练/验证/测试；预生成骨架划分尚未作为本轮训练实验 |
| `arrays.npz: baseline` | 5×978 | 给药前log2(TPM+1) |
| `drug_morgan` | 106×2048 | 0/1结构指纹 |
| `cell_index / drug_index` | 各3,364 | 每条样本取哪行背景/指纹 |
| `log10_dose_uM` | 3,364 | 以μM计的剂量取log10 |
| `expression_z24` | 3,364×978 | 24h实测响应答案 |
| `count72` | 3,364 | 72h处理/对照计数比答案 |
| `condition_weight` | 3,364 | 同条件跨库板权重合计为1 |
| `gene_id / library_condition_id` | 978 / 3,364 | 基因顺序和来源连接ID |
| `all_24h_candidates_with_QC.csv` | 3,410行 | 所有表达候选及46条排除记录的原因 |
| `instance_checks.csv / vehicle_plate_checks.csv` | 实例/板 | 原实例条件、缺失值和同板对照覆盖核对 |

## reference：前期已经核对的连接依据

这些表是由原实验整理出的中间依据，支持核查和重建当前训练集：

- `hms_conditions.csv`：4,140条库板条件的功能均值、背景/结构覆盖和表达签名对应。
- `hms_biorepeats.csv`：12,420条库板内生物重复，保留技术孔及源Excel行号。
- `lincs_signatures.csv / lincs_expression.npz`：当前输入覆盖内3,410条24h候选的说明和978基因数值，含之后排除的单实例条件。
- `baseline.csv`：按同一基因顺序排列的CCLE背景。
- `cells.csv / drugs.csv / input_features.npz`：固定的细胞、结构、基因和索引表；结构处理为原始SMILES解析、片段母体处理及含手性Morgan半径2指纹，化学分组字段另保留。
- `prior_condition_ids.csv`：用于核查早期配对准备条件与本轮仅需24h条件的差别，只有ID，没有重复的旧训练矩阵。

`training/audit.json`记录重建时实际读取的源文件路径和哈希；文件迁移、数组是否改变由`results/provenance.json`与根目录`MANIFEST.csv`说明。

## raw与sources：原文件和出处

`raw/`只包含三批HMS原Excel和LINCS小型说明表。大矩阵通过`download_manifest.csv/json`及下载脚本获取。`sources/`保留各来源数量与字段的核查依据；其中记录的本地规模和早期下载快照，不表示克隆后已经取得所有原矩阵。

原始孔、汇总签名、单细胞、候选连接和训练观测分别计数。其他来源当前用于后续设计与覆盖调查，不和本轮已训练样本混报。
