# HMS v3实验准备包

后续进展：此包原文件保持不变；已用于[第一部分的1次QC前联合模型调通训练/验证](../../../第一部分_依据整体框架的设计与起步_20261106/联合模型调通试验/run_summary.json)。外层测试未评分，生物学QC未冻结。下面“未运行模型拟合”是此包形成时的历史状态。

日期：2026-11-06。输入为已核对库板来源的v2，原文件和v2保持不变。此包完成矩阵组织、数值/索引/分组核查、描述性重复检查与划分草案；生物学QC与最终训练发布待后续完成。当前未运行模型拟合。

## 文件与数组

| 文件 | 内容 |
|---|---|
| arrays.npz | 按samples顺序的输入、目标、ID和权重，`allow_pickle=False`可读取 |
| samples.csv | 3407条库板条件的身份、剂量、来源、矩阵索引、原签名ID和条件分组 |
| cells.csv / drugs.csv / genes.csv | 五种细胞、106药、978基因及其数组索引/注释 |
| functional_replicates.csv | 10221条库板内生物重复、技术孔数量和原孔引用 |
| expression_provenance.csv | 6814条签名的组成实例、质量字段、原矩阵行和库板/条件 |
| split_assignments_DRAFT.csv | 四种泛化协议各五折＋两种剂量协议各一划分；每条记录的角色 |
| split_counts.csv | 各协议/折/角色的库板记录数和不同条件数 |
| cross_library_descriptive_checks.csv | 239个跨库板重复条件的功能差与表达相关 |
| audit.json | 数组、分组、化学、重复、缺失指标及软件版本的实际核查结果 |
| input_fingerprints.json | v2输入文件字节和SHA256，便于追踪来源变化 |

核心数组：

| 键 | 形状 | 含义 |
|---|---|---|
| expression_z | [3407,2,978] | L5响应z；第二维顺序3h/24h |
| expression_signature_id | [3407,2] | 对应实际签名ID |
| baseline_log2TPMp1 | [5,978] | CCLE未处理背景，基因顺序与响应相同 |
| drug_morgan | [106,2048] | 标准化母结构的radius=2、含立体Morgan位指纹 |
| cell_index / drug_index | [3407] | 每条观测对应的背景/结构索引 |
| log10_dose_uM | [3407] | log10浓度，原μM值保留在samples |
| relative_cell_count_72h | [3407] | 三个库板内生物重复均值的等权平均 |
| relative_net_growth_72h | [3407] | 同计数派生净增长，非独立测量 |
| count_bioreplicates | [3407,3] | 三个生物重复的技术孔均值，列文件名见bioreplicate_source |
| condition_weight | [3407] | 每条库板观测1/同条件库板数；总和3168 |
| expression_mask / functional_mask | [3407,2] / [3407] | 当前完整子集均有观测；不是通过生物学QC的标记 |

基因按Entrez ID而非旧符号对齐。保留原始SMILES、标准化含立体SMILES、保守划分用连接结构和Murcko骨架；用于划分的结构简化没有合并原化学身份或标签。

## 复现和验证

在项目根目录运行：

```bash
python3 data/scripts/package_hms_experiment_v3.py
python3 data/scripts/hms_ridge_baseline.py
python3 data/scripts/hms_ridge_baseline.py --task expression3
python3 data/scripts/hms_ridge_baseline.py --task expression24
```

后三条仅验证加载、索引、有限值与角色，默认不训练。拟合单任务单折的实验命令为：

```bash
python3 data/scripts/hms_ridge_baseline.py --fit --protocol unseen_drug --fold 0 --task function
```

上述`--fit`示例本轮未执行。执行后结果保存到独立baseline_pilots目录，标记为QC前试验；选择alpha只用验证集，输入标准化只拟合训练集。该脚本为独立岭回归，不是双任务共享表示模型，也不是完整论文实验。

本轮环境版本见audit.json，依赖numpy、pandas、scipy、scikit-learn、RDKit；从v2打包不需读取整个原始GCTX，不需GPU。

详见[项目总设计](../../../项目总设计_20261106/README.md)、[模型设计](../../../项目总设计_20261106/model_design.md)和[实验方案](../../../项目总设计_20261106/experiment_plan.md)。
