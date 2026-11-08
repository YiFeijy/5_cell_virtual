# 分子响应与功能读数的数据对接方法

这些方法分别回答不同问题，应在数据表中保留其证据等级。当前列出的12类方法覆盖本轮检索中主要可操作路线；它们不是可以随意互换的配对规则。

## 1 条件对应和实验配对不是一回事

| 等级 | 定义 | 能支持的分析 |
|---|---|---|
| A 协调设计的平行测量 | 同研究、确认细胞身份和药物、剂量与培养设计；表达和功能可在不同平行孔、不同终点测量 | 表达与功能之间的关系；在明确时间条件下预测后续终点 |
| B 跨平台精确条件对应 | 同细胞身份、同化学实体、同数值浓度；分别保留两个测量终点及assay | 跨实验关联与预测；不能按A等级解释 |
| C 浓度近似对应 | 在B的身份基础上允许有界剂量差，或曲线内插值 | 扩大关联样本；需要精确匹配子集和阈值敏感性分析 |
| D 只有细胞与药物对应 | 表达有具体剂量时间，功能只有IC50/AUC等曲线汇总 | 表达预测曲线敏感性；不支持浓度级活性真值 |
| E 不要求逐行对应 | 共享细胞或化学信息，但各任务分别有监督 | 联合训练或迁移学习；不能直接用未配对标签检验状态—功能关系 |

HMS本轮的3,852个条件是A等级候选：设计来自联合项目，但逐批次来源仍需复核。MIX-seq现成`1−AUC`主要属于D级功能标签。[HMS原始实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC5662764/)、[MIX-seq原始方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7453022/)。

## 2 十二类可用对接路线

| 方法 | 怎么做 | 文献或资源先例 | 使用条件与边界 |
|---|---|---|---|
| M1 同项目平行孔关联 | 追踪实验、细胞、药物、浓度、板和重复；分别记录表达时间与功能终点 | [HMS](https://doi.org/10.1038/s41467-017-01383-w) | 最优先；仍不是同个细胞的连续测量；一个后期功能标签不能因为多个表达时间而被当作多个独立标签 |
| M2 跨平台精确浓度关联 | 稳定细胞ID与完整化学身份映射后，匹配数值浓度，保留各自终点 | [LINCS与CTRP整合](https://doi.org/10.1093/nar/gkz805)中的严格子集 | 只解决身份和浓度，不自动消除培养、平台、批次和时间差异 |
| M3 最近浓度有界关联 | 同细胞药物，取最近log10浓度；例如差值小于0.2 | [Szalai代码](https://github.com/bence-szalai/Cell-death-signatures)与[WRFEN](https://doi.org/10.1186/s12859-020-03949-w) | 0.2约为1.58倍；建议同时报告精确、0.05、0.1、0.2阈值结果；不是所有药物都能容忍相同比例的剂量差 |
| M4 曲线内插值 | 在功能实测浓度范围内拟合或插值，得到表达实验剂量处的功能估计 | [scTherapy](https://doi.org/10.1038/s41467-024-52980-5) | 标签明确标为插值或拟合；报告拟合误差；无两侧支持的外推不进入主分析 |
| M5 表达到曲线敏感性 | 将特定时间的表达响应与同细胞药物的AUC/IC50关联 | [MIX-seq](https://doi.org/10.1038/s41467-020-17440-w)、[Dr.VAE](https://doi.org/10.1093/bioinformatics/btz158) | 功能任务是曲线汇总；不能把AUC复制成每个浓度的活性监督；曲线范围和方向要统一 |
| M6 双任务不完整标签训练 | 表达数据和功能数据分别计算对应损失，共享输入或部分表征 | [MultiDCP](https://doi.org/10.1371/journal.pcbi.1010367) | 不需要强行制造完整配对；联合收益须在真正有对应关系的独立队列上另验证 |
| M7 先学响应再迁移功能 | 先用扰动表达监督学习，再训练功能预测器 | [CSG²A](https://pmc.ncbi.nlm.nih.gov/articles/PMC11256952/) | 迁移设定不是数据配对；固定剂量或时间的模型输入不得写成真实缺失字段的恢复 |
| M8 预测响应作为功能特征 | 用预测表达、通路或潜在状态预测功能，与只用背景和药物比较 | [MultiDCP](https://doi.org/10.1371/journal.pcbi.1010367)、[TranSiGen](https://doi.org/10.1038/s41467-024-49620-3)提供相关表征学习路线 | 本项目需要训练外预测特征；不能让表达模型先见到功能测试集的对应表达标签 |
| M9 单细胞到条件级伪bulk | 按细胞系或细胞类型、药物、剂量、时间、实验重复聚合；与相同层级的功能标签连接 | [MIX-seq](https://doi.org/10.1038/s41467-020-17440-w)、[Chem-PerturBridge](https://github.com/theislab/Chem-PerturBridge) | 单细胞数不是独立功能样本数；保留每条件细胞数；不同生物学重复不要提前全部合并 |
| M10 统计分离功能相关表达 | 逐基因或通路回归功能读数，同时纳入背景及实验协变量，分析相关项和剩余项 | [MIX-seq](https://doi.org/10.1038/s41467-020-17440-w)、[Szalai](https://doi.org/10.1093/nar/gkz805) | 回归残差与潜在分支不自动具有因果可识别性；低表达质量不能解释成无响应 |
| M11 跨时间条件或分布关联 | 相同细胞背景、药物、剂量下连接早晚时间点；单细胞用群体分布或伪bulk | [NCI TPW](https://pmc.ncbi.nlm.nih.gov/articles/PMC6295263/)、[MIX-seq](https://doi.org/10.1038/s41467-020-17440-w)、[SequenTx](https://doi.org/10.1038/s42256-026-01192-1) | 破坏性测序不能提供同一细胞的真实轨迹；需要同批对照；单药时间外推与换药后的序贯状态转移分开验证 |
| M12 用辅助状态或功能数据迁移 | 以遗传扰动、形态、蛋白或肝毒性数据辅助表征；保留各自终点 | [ToxicoDB](https://academic.oup.com/nar/article/48/W1/W455/5838860)、[Way](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246468/)、[Szalai](https://doi.org/10.1093/nar/gkz805) | 共享靶点或通路只是关联线索；不把KO当成某个药物剂量，也不把形态变化当作直接死亡标签 |

M8中训练外特征和M3中的阈值比较是本项目的验证要求，不表示上述论文均执行了完全相同的设计。

## 3 建议的数据表

分开保存实验观测、稳定身份和跨观测关系，避免一行中把两个assay压成一个“样本”。

**表达观测表**

`expression_obs_id, source_release, experiment_id, plate_id, well_id, biological_replicate, cell_id, compound_entity_id, compound_batch_id, expression_dose_uM, expression_time_h, control_obs_ids, expression_representation, matrix_locator, quality_fields`

**功能观测表**

`functional_obs_id, source_release, experiment_id, plate_id, well_id, biological_replicate, cell_id, compound_entity_id, compound_batch_id, functional_dose_uM, functional_time_h, assay, readout_name, value, value_unit, normalization, label_origin, curve_id`

其中`label_origin`明确使用`measured`、`replicate_aggregate`、`interpolated`、`fitted`、`expression_derived`之一。缺失浓度的IC50/AUC汇总不填造假的`functional_dose_uM`。

**关系表**

`expression_obs_id, functional_obs_id_or_curve_id, pairing_grade, cell_mapping_evidence, compound_mapping_evidence, dose_log10_difference, expression_time_h, functional_time_h, same_study, same_experiment_design, ambiguity_flag`

## 4 身份映射的执行顺序

1. 细胞先保留原名，再对接Cellosaurus、DepMap、Sanger/COSMIC。删空格或连字符只生成候选，不能覆盖亲本、耐药株和基因改造株差异。[Cellosaurus官方说明](https://www.cellosaurus.org/faq)。
2. 药物优先使用来源提供的Broad ID或结构标识；保存完整InChIKey、原始SMILES及标准化结构。盐处理、立体异构体合并和Broad批次号折叠均记录依据。
3. 浓度先确认单位和数值字段。不要用LINCS舍入显示剂量覆盖更精确的实验剂量；DMSO百分数、蛋白质量浓度与小分子摩尔浓度不混算。
4. 时间分别保留两个终点。72小时不能填充到未知功能时间，也不能把24小时表达误标为72小时表达。
5. 基因用稳定ID映射，保留一对多关系和平台实际测量范围。bulk与单细胞伪bulk只在各自正确预处理后比较；978个landmark和推断基因分开。
6. 最后才生成关系表，报告精确、近似、歧义、缺失和未匹配数量；按细胞、药物、机制、剂量和时间检查交集偏倚。

## 5 划分和评价必须跟随数据层级

- 同一曲线、同条件技术重复、同孔细胞和派生伪bulk整体归入同一划分。
- 若评估未见细胞或药物，表达任务和功能任务都执行相应留出；同身份经过另一个任务进入训练也要说明。
- 剂量插值和时间插值与细胞、药物留出分别报告，不能用随机划分的分数代替。
- 先保留各数据集原始功能值和归一化方式；跨平台归一化不能让测试标签决定训练变换。
- 用实测早期表达的实验与用模型预测早期表达的实验分开。前者证明信息价值，后者才检验可部署的预测链。
