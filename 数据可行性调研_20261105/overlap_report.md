# 实际匹配规模与可执行任务

统计日期：2026年11月5日。所有数字由本地脚本重算。候选元数据对应、矩阵可读和可训练样本是不同阶段。

## Broad HMS 与 LINCS

| 项目 | 数量 | 解释 |
|---|---:|---|
| 功能原始表行数 | 12,420 | 含技术重复与重复条件 |
| 细胞系 | 6 | 含非肿瘤MCF10A |
| 药物Broad ID | 107 | 不代表全部已有结构输入 |
| 独立细胞药物组合 | 642 | 与剂量条件数区分 |
| 唯一功能剂量条件 | 3,852 | 全部在LJP表达元数据中找到条件对应 |
| 表达签名对应记录 | 9,355 | 跨来源与时间，不能当独立功能样本 |
| 同版本3小时和24小时条件 | 3,835 | GSE70138 |
| 同版本6小时和24小时条件 | 502 | GSE92742；不与上行直接相加 |
| 当前已取得结构字符串的药物 | 93 | 缺14个；字符串还未做化学有效性检查 |
| 有24Q4默认RNA profile索引的细胞 | 5 | 六细胞均有模型ID；MCF10A无该索引 |
| 同时具备上述结构与RNA索引的功能条件 | 2,790 | 还不是矩阵QC后的可训练数 |
| 同时具备上述输入及3/24小时记录的条件 | 2,781 | 结构补齐或更换基线来源后可能改变 |

表达端矩阵、控制与批次来源未完成核验，因此不报告完整联合预测训练样本数。

## MIX seq 原始表达试读

| 项目 | 数量或内容 |
|---|---|
| 药物 | Trametinib |
| 表达剂量 | 0.1 μM；来自原始补充表 |
| 时间 | 6小时和24小时，分别有DMSO对照 |
| normal单细胞合计 | 9,894 |
| 四个条件各至少20个normal细胞的细胞系 | 24 |
| 上述细胞系且有作者敏感性标签 | 24 |
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
