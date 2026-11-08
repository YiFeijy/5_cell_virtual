# MCF10A在当前CCLE矩阵中缺行的原因追查

日期：2026年11月6日。

## 结论分成三个层次

**已经证实：**MCF10A的ModelID没有找错，本地文件与官方文件一致；DepMap 24Q4的模型名录登记了它，但该发布的omics profile索引和默认profile索引没有为它提供RNA记录。因此，当前选用的RNA矩阵没有它，是这套固定发布的RNA覆盖情况，而不是查错拼写或本地文件截断。

**仍未查明：**发布说明没有说明这个具体模型为何没有RNA profile。仅凭这些公开文件，无法区分“该项目未产生相应RNA测量”“测过但未进入该公开发布”等上游情形。这里不补写一个未经证实的实验原因。

**已补充的事实：**MCF10A确实有其他表达数据。本地LJP矩阵中已读取295个MCF10A的DMSO对照表达；另已下载并验证GSE152410中的四个0小时MCF10A RNA-seq对照样本。当前缺的是所选CCLE版本内的RNA行，不是整个研究领域没有这个细胞的表达数据。

我此前在发现缺行后没有继续把原因与替代来源查清，说明停在了检索结果层面。这里补齐证据链，并将对异常追问原因写入项目必读要求。

## 1 是否下载损坏或读错文件

重新遍历本地文件计算MD5，与官方24Q4文件清单逐一比较：

| 文件 | 检查结果 |
|---|---|
| Model.csv | 官方MD5一致 |
| OmicsProfiles.csv | 官方MD5一致 |
| OmicsDefaultModelProfiles.csv | 官方MD5一致 |
| OmicsExpressionProteinCodingGenesTPMLogp1.csv | 官方MD5一致 |
| 本轮补充的README、ModelCondition、默认条件profile表 | 官方MD5一致 |

这证明这些本地文件与固定发布文件的字节内容一致。数据完整性核查与生物学质量评价分别记录。详细校验值见[文件完整性复核](CCLE_file_integrity_recheck.json)。

## 2 是否细胞名字或模型ID找错

在Model.csv中按原名、去标点名、CCLE名和稳定ID查找，得到同一条模型记录：

| 字段 | 值 |
|---|---|
| CellLineName | MCF 10A |
| StrippedCellLineName | MCF10A |
| CCLEName | MCF10A_BREAST |
| ModelID | ACH-001357 |
| Cellosaurus RRID | CVCL_0598 |

再以ACH-001357查询表达矩阵的行ID，匹配数为0。查找范围覆盖了名称与稳定ID。[实际模型记录](MCF10A_DepMap_model_record.csv)

## 3 是否只是选了蛋白编码基因表达表导致缺行

本轮继续检查了相同发布中的映射资料：

| 检查对象 | ACH-001357相关记录数 | 说明 |
|---|---:|---|
| Model.csv | 1 | 登记了模型身份 |
| OmicsProfiles.csv | 0 | 全部公开omics profile映射中没有该模型条目 |
| OmicsDefaultModelProfiles.csv | 0 | 没有默认模型profile |
| ModelCondition.csv | 0 | 没有该模型条件条目 |
| 所选RNA表达矩阵 | 0 | 没有该模型表达行 |

[官方README](DepMap24Q4_README.txt)说明，Model.csv登记的是DepMap门户中各类数据引用的模型；OmicsProfiles用于profile与模型的映射；默认profile表指出模型级数据选择了哪个profile。模型名录与RNA数据覆盖服务不同范围。

本地Model.csv有2,105个模型，表达矩阵有1,673行。一个模型在名录中有身份记录，并不等于这个发布为它提供了RNA测量。对MCF10A而言，缺失已定位到profile映射这一层，而不仅是某个基因列表的选择。

## 4 是否因为它是非癌细胞

仅靠这个解释不成立。官方Model.csv标为Non-Cancerous的137个模型中，85个实际出现在这张RNA表达矩阵中。[分类与表达覆盖表](NonCancerous_model_expression_coverage.csv)

因此，细胞类别本身不足以解释MCF10A为何未纳入该发布。具体的实验或发布决定仍须有来源证据；本轮已下载的README未给出该模型的单独说明。

## 5 现有本地数据中是否已有MCF10A参照表达

有。本地LINCS Phase II的LJP实例说明中有295个MCF10A的vehicle对照：3h为149个，24h为146个，来自12个检测板。

这295个实例ID均在Level3矩阵中找到，实际提取为295×978的实测landmark表达；读取值均有限，未出现`-666`缺失标记。[实例说明](MCF10A_LJP_vehicle_instances.csv)、[实际表达矩阵](MCF10A_LJP_vehicle_landmark_values.npz)

这些值可以支持“同平台同期参照状态”这一输入定义。其时间为3h或24h，命名与时间字段保留原义。若使用这种参照，可为所有六种细胞统一建立同平台输入，再与CCLE背景输入的方案比较。

## 6 是否存在公开RNA-seq的零时点样本

已从GEO GSE152410下载四个实际文件：

| GEO样本 | 样本标题 | 实验注释 |
|---|---|---|
| GSM4614840 | ctrl_0_C1_A_RNAseq | MCF10A、PBS、0 ng/ml、0h |
| GSM4614841 | ctrl_0_C1_B_RNAseq | MCF10A、PBS、0 ng/ml、0h |
| GSM4614842 | ctrl_0_C1_C_RNAseq | MCF10A、PBS、0 ng/ml、0h |
| GSM4614843 | ctrl_0_C1_D_RNAseq | MCF10A、PBS、0 ng/ml、0h |

四个文件均为Kallisto转录本丰度表，包含target_id、length、eff_length、est_counts、TPM。已按GEO文件清单核对大小，并完整读取gzip检查CRC；每个文件199,169条转录本记录，TPM均有限且非负。

根据target_id自带的Ensembl基因标识，将转录本TPM按基因求和，保留四个样本各自的数值，得到60,554行基因注释。[四个样本的基因TPM](MCF10A_RNAseq_GSE152410/MCF10A_gene_TPM_four_controls.csv)、[样本和完整性记录](MCF10A_RNAseq_GSE152410/file_content_checks.json)

这批数据来自MCF10A多组学配体研究。作者方法显示，实验培养基去除了EGF和insulin，与常规增殖培养有差别。因此它是需要核对培养背景和基因ID的独立RNA-seq候选来源，当前没有被直接插进CCLE矩阵。[原始研究及数据说明](https://pmc.ncbi.nlm.nih.gov/articles/PMC9546880/)、[GEO归档](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152410)

按原始基因符号进行一对一候选匹配时，974个landmark有唯一候选；PIK3R3和PAK6存在同符号多Ensembl条目，HDGFRP3与PLSCR3未按原符号得到唯一候选。这些是ID/注释复核项，尚未用零值补齐或作为最终978维输入。[候选基因对应](MCF10A_RNAseq_GSE152410/landmark_gene_symbol_mapping_candidates.csv)

## 7 这对当前对接规模意味着什么

之前的3,168个条件，是“固定使用当前CCLE版本、结构字符串可解析、表达和功能时点齐全”这一输入选择下的候选规模。这个数字不是MCF10A可用性的生物学结论，也不是全部数据的最终规模。

MCF10A目前有637个同时找到3h和24h表达的功能条件。其记录继续保留。接下来的工作应比较两种清楚定义的输入：统一来源的细胞背景RNA，以及所有六细胞的同平台同期vehicle参照。新找到的RNA-seq需要先核对培养与注释，采用哪种输入及如何校准，要以实际比较结果决定。

## 8 证据与复现

- [机器可读原因追查与进度](MCF10A原因核查与对接进度.json)。
- [复现脚本](../scripts/investigate_mcf10a_and_progress.py)。
- [官方24Q4发布](https://doi.org/10.25452/figshare.plus.27993248.v1)。
- 原始发布文件未修改；本轮补充的独立来源数据保存在本调查目录。
