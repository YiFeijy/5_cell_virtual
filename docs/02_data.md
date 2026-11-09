# 02 数据是什么、怎样测量、为什么使用

更新：2026-11-09。先理解每种数据的实验和指标，再讨论怎么连接。本页把**研究介绍、数据发布、文件直链**分别列出；细胞覆盖来自本地实际观测表/矩阵，不用一张通用细胞名单替代。

## 先分清本项目里的四种数量

1. **实验/发布范围**：论文或官网介绍的总体设计。
2. **本地实际覆盖**：下载版本中真正出现的模型、样本、孔、细胞或签名。
3. **可核查的连接关系**：按身份、剂量等找到的待复核对应，不新增实验。
4. **本轮训练范围**：CCLE＋LINCS Phase II＋HMS形成的5细胞、106药、3364条库板观测/3133个不同条件。

下面“细胞数”若指单细胞记录，会明确写细胞/核；“细胞系/模型数”则按源ID统计。全量细胞名单位于[data/catalog](../data/catalog/README.md)，保留命名和ID，不凭相似拼写合并。

## 按来源阅读

| 数据 | 本地实际细胞覆盖 | 测量内容 | 当前作用 |
|---|---|---|---|
| [CCLE / DepMap 24Q4](#ccle_24q4) | 1673个有RNA的ModelID | 背景RNA-seq | 本轮输入 |
| [LINCS Phase II / GSE70138](#lincs_gse70138) | 41个cell_id；30个有小分子签名 | L1000响应 | 本轮表达答案 |
| [HMS LINCS Joint Project](#hms) | 6种，完整名单在正文 | 成像核计数 | 本轮功能答案 |
| [CTRPv2](#ctrpv2_2015) | 2015实测887；2016 AUC 835、表达/CNV各823 | ATP相关活性 | 功能扩展 |
| [GDSC1](#gdsc1) | 970个源名称/模型ID | 染料/代谢信号 | 功能扩展 |
| [GDSC2](#gdsc2) | 969个源名称/模型ID | ATP相关强度 | 功能扩展 |
| [PRISM Repurposing 19Q4](#prism) | 初筛578行/568个有当前名称；复筛489行/481个有当前名称 | 相对条形码丰度 | 功能扩展 |
| [LINCS Phase I / GSE92742](#lincs_gse92742) | 76个cell_id；71个有小分子签名 | L1000响应 | 额外表达监督 |
| [MIX-seq](#mixseq) | 169个normal且有DepMap ID的背景 | 单细胞RNA | 表达外部验证和异质性 |
| [sci-Plex2](#sciplex2) | A549 | 单细胞/核RNA | 小规模表达流程验证 |
| [sci-Plex3](#sciplex3) | A549、K562、MCF7，另有缺注释行 | 单细胞/核RNA | 表达外部验证与组成研究 |
| [NCI TPW / GSE116436](#nci_tpw) | NCI样本中实际60个来源 | 芯片表达 | 时间与表达外部验证 |
| [Novartis DRUG-seq U2OS MoABox](#drugseq) | U2OS | 每孔UMI计数 | 扩大药物表达覆盖 |
| [Tahoe-100M](#tahoe_shard) | 表达分片50；注释表102个不同身份 | 单细胞RNA counts | 格式核对与选择性扩展 |

<a id="ccle_24q4"></a>

## 1. CCLE / DepMap 24Q4：细胞基线RNA

**介绍在哪里，谁产生的数据**

Broad Institute的CCLE/DepMap用于系统刻画细胞模型的分子背景。当前使用的是24Q4 Public固定发布，不是将CCLE早期论文中的模型数直接当成当前文件规模。


[项目介绍](https://depmap.org/portal/data_page/) · [固定24Q4发布](https://doi.org/10.25452/figshare.plus.27993248.v1) · [CCLE方法论文](https://doi.org/10.1038/s41586-019-1186-3)


**有哪些细胞系，实际有多少数据**

实际表达文件有1,673个ModelID、19,193个蛋白编码基因列；模型说明表有2,105条模型记录，二者覆盖不同。细胞名称和组织直接按ModelID从Model.csv连接，完整名单列出每个RNA模型，而不是只列肿瘤类型。

本轮用到SKBR3、MCF7、HS578T、BT20、MDAMB231；同文件还包含A549、PC3、HT29等多种背景，完整名单见下方CSV。MCF10A有模型身份记录，但固定24Q4没有相应RNA profile，原因追查见[data_identity.md](data_identity.md)。

从实际覆盖中举例：MCF7、MDA-MB-231、SK-BR-3、A549、PC-3、K-562、1156QE/8、1273/99、143B、1618-K、170-MG-BA、1777N Rpmet。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/ccle_24q4.csv)。


**实验怎样测出来**

从常规培养细胞取得RNA，经RNA测序、基因表达定量及TPM处理后发布。TPM描述按转录本长度和样本总量校正的相对丰度；本文件再存为log2(TPM+1)。它测的是独立实验中的背景RNA，不是HMS每个给药孔的起始样本。


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| ModelID | 稳定模型身份，如ACH-000536对应BT20 |
| log2(TPM+1) | 表达相对丰度的对数变换；0对应TPM=0 |
| 基因列“符号 (Entrez ID)” | 定位基因，跨来源对齐使用ID |
| ProfileID / ProfileType | 表明模型具体有哪些组学测量；Model表本身不保证有RNA |

这里没有药物剂量或给药时长。背景RNA可作为输入，其测量定义与LINCS响应z值分别保留。


**下载了什么，从哪里下载**

表达矩阵的一行是一个细胞模型、一列是一个基因。Model.csv解释身份和组织；OmicsDefaultModelProfiles.csv及OmicsProfiles.csv解释模型与RNA等测量的对应。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | [下载文件](https://ndownloader.figshare.com/files/51065489) |
| `Model.csv` | [下载文件](https://ndownloader.figshare.com/files/51065297) |
| `OmicsDefaultModelProfiles.csv` | [下载文件](https://ndownloader.figshare.com/files/51065339) |
| `OmicsProfiles.csv` | [下载文件](https://ndownloader.figshare.com/files/51065723) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

本轮提取五细胞与L1000共同的978个实测基因，提供给药前背景。后续扩细胞时先检查真实RNA覆盖、身份和培养条件，再加入输入。


[回到来源导航](#按来源阅读)


<a id="lincs_gse70138"></a>

## 2. LINCS Phase II / GSE70138：药物扰动表达

**介绍在哪里，谁产生的数据**

NIH LINCS计划中Broad的L1000/CMap资源，目的是记录细胞对化合物或遗传扰动的转录响应，比较作用模式及分子联系。


[研究与数据发布页](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) · [实验和层级说明](https://clue.io/connectopedia/data_levels) · [原始研究](https://doi.org/10.1016/j.cell.2017.10.049)


**有哪些细胞系，实际有多少数据**

sig_info实际出现41个原始cell_id，其中30个具有小分子trt_cp签名。这里还含遗传扰动、对照、原代/分化模型及带后缀的模型；不把41个ID全部当成同一类癌细胞系，也不把Touchstone常用九细胞当成本文件全部覆盖。

HMS对应的六个是BT20、HS578T、MCF10A、MCF7、MDAMB231、SKBR3。带.101、.311等后缀的cell_id按发布原样保留，先查cell_info中的modification和base_cell_id再决定关系。

从实际覆盖中举例：BT20、MCF7、HS578T、MDAMB231、SKBR3、A549、PC3、HT29、A375、ASC、ASC.C、CD34。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/lincs_gse70138.csv)。


**实验怎样测出来**

培养细胞按药物、剂量和时间处理后提取RNA。L1000使用靶向探针、扩增及Luminex多重荧光读出测量978个landmark转录本；其他特征由模型推断。发布流程将原始信号处理成归一化表达，再转成标准化差异，并汇总重复形成签名。方法见[原论文](https://doi.org/10.1016/j.cell.2017.10.049)及[CLUE说明](https://clue.io/connectopedia/data_levels)。


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| Level3 | 单个实例的归一化表达，用于本轮检查实例数值、身份与对照覆盖 |
| Level5 / MODZ | 重复汇总的标准化响应；正负表示相对参考的上调/下调，不是TPM或直接倍数 |
| sig_id / distil_id | 汇总签名及其原始实例ID链 |
| pert_id / pert_idose / pert_itime | 扰动身份、带单位剂量和处理时间 |
| pert_type | 区分小分子trt_cp、vehicle和其他实验 |
| TAS / distil_cc_q75 | 转录信号及重复一致性相关质量信息；不是活性或死亡指标 |

本轮固定24h及实际0.04—10 μM网格。全发布时间/剂量依各行记录；溶剂对照和检测板用inst_info定位。本轮不按TAS强弱排除。


**下载了什么，从哪里下载**

Level3为345,976个实例×12,328特征，Level5为118,050条签名×12,328特征；107,404条属于小分子trt_cp签名。sig_info说明条件与重复，inst_info追实例/板，gene_info标记实测基因，pert_info解释药物及结构。发布总量与本轮3364条观测分别使用。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70138&format=file&file=GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz) |
| `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70138&format=file&file=GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz) |
| `GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70138&format=file&file=GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz) |
| `GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70138&format=file&file=GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

按库板、细胞、药物和实际剂量与HMS建立对应，本轮提取978维24h响应作答案。对接步骤见[03](03_integration.md)。其余表达可用于额外监督，但先核对重复与测试药物隔离。


[回到来源导航](#按来源阅读)


<a id="hms"></a>

## 3. HMS LINCS Joint Project：原孔成像计数

**介绍在哪里，谁产生的数据**

HMS LINCS与Broad联合项目，Hafner、Niepel、Sorger等团队在相应细胞/药物条件下平行测量转录与生长表型。数据集20237/20238/20239分别提供一批实验。


[项目数据集20237](https://lincs.hms.harvard.edu/db/datasets/20237/) · [第二批20238](https://lincs.hms.harvard.edu/db/datasets/20238/) · [第三批20239](https://lincs.hms.harvard.edu/db/datasets/20239/) · [原始研究](https://doi.org/10.1038/s41467-017-01383-w)


**有哪些细胞系，实际有多少数据**

六种细胞完整名单是：**BT20、HS578T、MCF10A、MCF7、MDAMB231、SKBR3**。MCF10A为非肿瘤乳腺上皮模型，其余五种为乳腺癌来源。本轮依据固定CCLE背景覆盖，采用除MCF10A外的上述五种；MCF10A的计数和L1000观测仍保留，背景来源另查。

三个文件共36,912条孔记录；107个药物×6剂量×6细胞对应3,852个不同条件。保留LJP005/006库板时有4,140条库板条件，再保留三批生物重复为12,420条重复记录。

从实际覆盖中举例：BT20、HS578T、MCF10A、MCF7、MDAMB231、SKBR3。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/hms.csv)。


**实验怎样测出来**

原Excel的Meta表描述：铺板培养24h后给药；平行板测给药开始的细胞数量；药物暴露72h后固定，以Hoechst染色和图像扫描识别细胞核，估计每孔细胞数。DMSO孔作为终点参照。技术孔先在每批、每库板条件内平均，再等权汇总三批。原实验说明在各Excel的Meta工作表，逐字段解释在Data Columns工作表，原孔记录在Data工作表。


**指标和字段怎样理解**

| 原字段 | 定义和用途 |
|---|---|
| Total Cell Count Before Treatment | N0，给药开始的平行孔计数 |
| Total Cell Count After Treatment | Nt，72h处理组计数 |
| Total Control Cell Count | Nc，72h DMSO参照计数 |
| Relative Cell Count | Nt/Nc；本轮功能答案，1表示与对照相当，0.75表示数量为对照75% |
| Relative Growth Rate After Treatment | 原表数值对应(Nt−N0)/(Nc−N0)，描述相对净增长；公式与Hafner归一化GR指标不同 |

剂量为0.04、0.12、0.37、1.11、3.33、10 μM，终点为72h。计数比包含生长减慢和数量减少的共同后果，不直接改称死亡细胞比例。


**下载了什么，从哪里下载**

三份Screen20237/20238/20239 Excel都有实验说明及Data表，每行是一个孔。Data中保留细胞、BRD药物ID、库板、孔位、重复、浓度/单位和计数。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `Screen20237_LJP_processedData1.xls` | [下载文件](https://lincs.hms.harvard.edu/_static/db/prod-20200624/Screen20237_LJP_processedData1.xls) |
| `Screen20238_LJP_processedData2.xls` | [下载文件](https://lincs.hms.harvard.edu/_static/db/prod-20200624/Screen20238_LJP_processedData2.xls) |
| `Screen20239_LJP_processedData3.xls` | [下载文件](https://lincs.hms.harvard.edu/_static/db/prod-20200624/Screen20239_LJP_processedData3.xls) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

它与GSE70138有项目库板及实验条件可核对，因此先用于来源清楚的联合比较。72h计数和24h表达是平行实验的不同终点，保留其关系而不虚构单细胞轨迹。例子与重算见[03](03_integration.md)和src/aggregate_hms.py。


[回到来源导航](#按来源阅读)


<a id="ctrpv2_2015"></a>

## 4. CTRPv2：浓度级ATP活性与药敏曲线

**介绍在哪里，谁产生的数据**

Broad/CTD²的Cancer Therapeutics Response Portal用于研究癌细胞的分子特征与化学敏感性，寻找癌细胞对特定药物的依赖。


[项目介绍](https://portals.broadinstitute.org/ctrp.v2/) · [2015原研究](https://doi.org/10.1158/2159-8290.CD-15-0235) · [原包URL的维护配置](https://github.com/BHKLAB-DataProcessing/CTRPv2-Pharmacoset_Snakemake/blob/main/config/pipeline.yaml) · [整理版发布](https://zenodo.org/records/4429239)


**有哪些细胞系，实际有多少数据**

本地2015扩展表实际涉及887个细胞ID，完整名单由浓度表的experiment_id连接实验表再连接细胞表得到。细胞说明表有1,107行，不能直接用它代替实际药物测量覆盖。

本地有6,171,005条浓度记录、545个处理条目；其中49个明确为固定比例组合，496个未带该组合标记。官网核心筛选介绍的481化合物/860细胞，与这个扩展文件的覆盖分别标注。

从实际覆盖中举例：BT20、MCF7、HS578T、MDAMB231、SKBR3、A549、HCT116、PC3、K562、HT29、U2OS、VCAP。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/ctrpv2_2015.csv)。

2016补充包另外核查：AUC表有326,148条曲线记录、835个细胞ID；基线表达和拷贝数表各覆盖823个细胞。完整名单保留每个ID在三类表中的覆盖：[2016细胞清单](../data/catalog/ctrpv2_2016.csv)。


**实验怎样测出来**

细胞在多个药物浓度下培养，72h后用CellTiter-Glo读出ATP相关发光信号，经过背景和对照处理及重复汇总，形成相对活性，并拟合剂量曲线。原研究使用16点浓度系列；各记录实际浓度仍按cpd_conc_umol保存。方法来源是[Seashore-Ludlow等原研究](https://doi.org/10.1158/2159-8290.CD-15-0235)及下载包README/COLUMNS。


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| cpd_conc_umol | 实际浓度，μM |
| cpd_avg_pv | 重复汇总的实测相对活性；1为参照量级，是比例而非百分数 |
| cpd_pred_pv | 拟合曲线在该浓度的值，与实测值分列 |
| AUC及曲线参数 | 概括一个细胞—药物实验的浓度系列，不是单一剂量的计数 |
| master_ccl_id / experiment_id / master_cpd_id | 细胞、实验和处理身份 |

当前实测相对值范围约0.001922—2.881；大于1的值保留。原孔、汇总浓度和曲线分别计数。

**2016补充的指标也要分开：**area_under_curve是拟合浓度曲线的积分，本地范围0.0691—29.35，不是通用0—1活性；mrna_expression_avg_log2为基线平均log2表达，copy_number_var_avg_log2为相对log2拷贝数。它们提供药敏与背景特征，表达不是给药后的响应。

**PharmacoSet已实际读取：**cell对象887行，其master_ccl_id集合与2015实测细胞集合相同；drug对象544行、sensitivity info为395,263条曲线。544与原545处理条目的差别已定位：整理对象将tipifarnib-P2/P1的源ID417416、417417并为一条417416///417417。后续继续保留两个原ID和结构核对，不能把整理行数变化当成新增或丢失实验。profiles内还分published/recomputed AAC/EC50等字段，应按对象字段定义使用。[本次对象核查](../data/catalog/pharmacoset.json)。


**下载了什么，从哪里下载**

2015 ZIP包含per_cpd_well原孔、post_qc浓度表、curves_post_qc及细胞/药物/实验说明。2016补充ZIP包含AUC、基线表达、拷贝数及关联资料；PharmacoSet RDS是同一实验来源的整理版，不作为一份新的独立实验重复。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `CTRPv2.0_2015_ctd2_ExpandedDataset.zip` | [下载文件](https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip) |
| `CTRPv2.1_2016_pub_NatChemBiol_12_109.zip` | [下载文件](https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.1_2016_pub_NatChemBiol_12_109.zip) |
| `CTRPv2_PharmacoSet.rds` | [下载文件](https://zenodo.org/api/records/4429239/files/CTRPv2.rds/content) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

用于扩展功能监督。先拆分单药与组合，按细胞ID与结构核对，再保留原浓度、时间和实验来源连接表达。已有146条“数值剂量相等”的初筛关系涉及100个表达查询，仅是待核对关系；本轮训练没有使用CTRP。


[回到来源导航](#按来源阅读)


<a id="gdsc1"></a>

## 5. GDSC1：染料/代谢信号及拟合药敏

**介绍在哪里，谁产生的数据**

Sanger与MGH联合产生的Genomics of Drug Sensitivity in Cancer第一代筛选，用于关联细胞分子特征与药敏。


[项目与测量方法](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/) · [数据发布页](https://cellmodelpassports.sanger.ac.uk/downloads/gdsc) · [原始研究](https://doi.org/10.1016/j.cell.2016.06.017)


**有哪些细胞系，实际有多少数据**

27Oct23原始表实际有970个不同细胞名称及970个非空Sanger模型ID；覆盖名单由原始孔表直接统计，不从另一版官网清单替代。共5,976,183条孔记录，含背景和对照；其中library单药孔3,179,341条、403个药物ID。

从实际覆盖中举例：MCF7、MDA-MB-231、A549、HCT-116、PC-3、K-562、22RV1、23132-87、42-MG-BA、451Lu、5637、639-V。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/gdsc1.csv)。


**实验怎样测出来**

在板中铺细胞并施加剂量系列，72h后以Syto60染料或Resazurin代谢信号读出。先结合板内对照处理强度，再拟合反应曲线。官方方法列出9点二倍或5点四倍稀释设计；具体样本按原孔浓度与ASSAY保存。[官方方法与字段说明](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/)


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| INTENSITY / ASSAY | 原始检测强度及检测类型；强度还未等同相对活性 |
| TAG / BARCODE / CONC | 孔类型、板ID、原药物浓度；归一化需区分处理和对照 |
| SANGER_MODEL_ID / COSMIC_ID | 保留源身份，跨库对接时核对映射 |
| LN_IC50 | 拟合IC50（以μM计）的自然对数；不是某孔的存活率 |
| AUC | 规定浓度范围内的拟合曲线面积比例 |
| RMSE | 曲线拟合误差，反映拟合与读数的一致程度 |

原孔读数与fitted文件各自保留。研究给药72h；原表部分DURATION缺失，按源记录与实验说明分别标注。


**下载了什么，从哪里下载**

raw ZIP解压后为孔级CSV及说明；fitted XLSX保存细胞—药物曲线结果。不同位点/批次可能有同一化合物的不同Drug ID，因此后续按结构和来源复核，而不只靠名称合并。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `GDSC1_public_raw_data_27Oct23.zip` | [下载文件](https://cmp.cog.sanger.ac.uk/download/GDSC1_public_raw_data_27Oct23.zip) |
| `GDSC1_fitted_dose_response_27Oct23.xlsx` | [下载文件](https://cmp.cog.sanger.ac.uk/download/GDSC1_fitted_dose_response_27Oct23.xlsx) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

用来扩大功能条件及检验检测方法变化的影响。当前已核原孔和对照，下一步执行来源内归一化及少量原实例复查，再以保留剂量的功能表接入部分标签学习。本轮未使用其标签。


[回到来源导航](#按来源阅读)


<a id="gdsc2"></a>

## 6. GDSC2：ATP相关信号及拟合药敏

**介绍在哪里，谁产生的数据**

Sanger采用更新培养板及筛选流程的第二代GDSC。与GDSC1同属药敏资源，但检测流程和批次分别保存。


[项目、方法和字段](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/) · [发布页](https://cellmodelpassports.sanger.ac.uk/downloads/gdsc)


**有哪些细胞系，实际有多少数据**

27Oct23原孔中有969个不同细胞名称、969个非空Sanger模型ID，15,069,449条孔记录。其中library单药孔3,166,905条、296个药物ID。官方当前说明页列297化合物，与本地这个版本及孔标签筛选下的药物ID数分别报告。

从实际覆盖中举例：MCF7、MDA-MB-231、A549、HCT-116、PC-3、K-562、22RV1、23132-87、42-MG-BA、451Lu、5637、639-V。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/gdsc2.csv)。


**实验怎样测出来**

官方流程为1536孔板铺板24h后给药，药物处理72h后以CellTiter-Glo测ATP相关发光。采用7点浓度系列，依板内DMSO和背景对照处理原强度后拟合。[测量和处理定义](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/)


**指标和字段怎样理解**

主要字段与GDSC1相同：INTENSITY是原强度，CONC是浓度，TAG是孔类型，LN_IC50/AUC是拟合结果。两代的实验来源和ASSAY保留。

**本地DURATION全为4，而官方描述给药72h。** 官方流程另有铺板后24h再给药；原导出字段与给药时长是不同层的记录，不能把4直接重命名为96h药物暴露。此处分别保留原字段与官方协议说明，后续按实验/板追查其具体计时口径。


**下载了什么，从哪里下载**

raw ZIP内有孔级CSV；fitted XLSX另有细胞、药物及曲线指标。全部孔数包含控制及未用孔，不是15,069,449个独立细胞—药物条件。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `GDSC2_public_raw_data_27Oct23.zip` | [下载文件](https://cmp.cog.sanger.ac.uk/download/GDSC2_public_raw_data_27Oct23.zip) |
| `GDSC2_fitted_dose_response_27Oct23.xlsx` | [下载文件](https://cmp.cog.sanger.ac.uk/download/GDSC2_fitted_dose_response_27Oct23.xlsx) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

适合扩展浓度级功能监督，并和GDSC1/HMS分别比较来源影响。已完成孔类型及板对照覆盖核查；后续先处理强度，再核对细胞身份、药物结构和剂量。当前没有将其强度直接作为HMS计数比。


[回到来源导航](#按来源阅读)


<a id="prism"></a>

## 7. PRISM Repurposing 19Q4：混合细胞的条形码丰度

**介绍在哪里，谁产生的数据**

Broad的PRISM Repurposing用于大范围药物重定位筛选，检测已有化合物对不同癌细胞的选择性作用。设计介绍为初筛4,518药、复筛1,448药；下载矩阵的处理列含剂量与筛选批次，列数不是唯一药物数。


[研究介绍、方法和下载入口](https://depmap.org/repurposing/) · [原始研究](https://doi.org/10.1038/s43018-019-0018-6)


**有哪些细胞系，实际有多少数据**

必须把矩阵行与认证细胞分开：初筛矩阵578行中，568行有当前细胞名称，10行带FAILED_STR；复筛489行中，481行有当前名称，8行带FAILED_STR。完整清单同时保存原矩阵行ID、当前DepMap ID、原分配名称、STR标记和原备注。

例如ACH-000028_FAILED_STR原先标为KPL1，原备注显示STR匹配MCF7；本次保留该原始异常行，不把它直接并为MCF7。初筛细胞说明表588行，也不是初筛矩阵实际578行。

[初筛完整行/细胞表](../data/catalog/prism_primary.csv)、[复筛完整行/细胞表](../data/catalog/prism_secondary.csv)；含FAILED_STR标记及原备注。


**实验怎样测出来**

先给不同细胞系引入可区分的稳定分子条形码，将多个细胞系混合培养后给药；五天后裂解，通过条形码捕获及Luminex读出估计各细胞系的相对丰度。信号相对vehicle处理并做发布流程的重复/批次处理。[原论文方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7328899/)及[官方README入口](https://depmap.org/repurposing/)。


**指标和字段怎样理解**

| 指标/字段 | 意义 |
|---|---|
| replicate-collapsed logfold-change | log2相对条形码丰度，0对应比值1、−1对应比值0.5；独立于RNA表达 |
| dose / broad_id / screen_id | 浓度、药物样品与筛选来源 |
| depmap_id / row_name | 细胞身份与实际矩阵行；FAILED_STR行另外保留 |
| fitted AUC / IC50 | 整条曲线的汇总，与浓度级实测相对丰度分开 |

初筛矩阵578×4,686；复筛489×13,008，其中5,874,229个非缺失数值。复筛设计是从10 μM开始的八点四倍稀释。五天读数与本轮72h计数分别标注。


**下载了什么，从哪里下载**

已下载两张重复汇总矩阵、各自细胞/处理说明、README和复筛曲线参数。原始MFI等官方另有发布，本批并未将所有中间文件都下载。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `primary-screen-replicate-collapsed-logfold-change.csv` | [下载文件](https://ndownloader.figshare.com/files/20237709) |
| `primary-screen-cell-line-info.csv` | [下载文件](https://ndownloader.figshare.com/files/20237718) |
| `secondary-screen-replicate-collapsed-logfold-change.csv` | [下载文件](https://ndownloader.figshare.com/files/20237757) |
| `secondary-screen-cell-line-info.csv` | [下载文件](https://ndownloader.figshare.com/files/20237769) |
| `secondary-screen-replicate-collapsed-treatment-info.csv` | [下载文件](https://ndownloader.figshare.com/files/20237763) |
| `secondary-screen-readme.txt` | [下载文件](https://ndownloader.figshare.com/files/20238123) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

可做不同药物范围及不同功能实验的外部比较。按经过STR审核的身份、药物样品和实际剂量连接，FAILED_STR记录单独调查，不去掉后缀后冒充认证模型。当前未进入训练。


[回到来源导航](#按来源阅读)


<a id="lincs_gse92742"></a>

## 8. LINCS Phase I / GSE92742：较早发布的扰动表达

**介绍在哪里，谁产生的数据**

同一Broad L1000资源的Phase I固定发布。研究目的和测量平台与Phase II相通，但实验范围、批次和模型覆盖不同，后续按具体记录关联。


[发布和研究说明](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) · [平台方法](https://doi.org/10.1016/j.cell.2017.10.049)


**有哪些细胞系，实际有多少数据**

实际签名说明表覆盖76个原始细胞模型ID，其中71个有小分子trt_cp签名。它包含A375、A549、BT20、MCF7、MCF10A、MDAMB231、SKBR3、PC3、VCAP、U2OS等，也含原代/分化模型。完整清单保留每个模型的类型和组织注释。

从实际覆盖中举例：BT20、MCF7、HS578T、MDAMB231、SKBR3、A549、HCT116、PC3、HT29、VCAP、A375、A673。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/lincs_gse92742.csv)。


**实验怎样测出来**

同样以L1000的靶向探针和多重荧光读出测量landmark表达，随后做归一化、标准化响应及重复汇总。具体时间、剂量和扰动类型按其sig_info/inst_info读取，不从Phase II借用统一设置。平台方法见[Cell原始论文](https://doi.org/10.1016/j.cell.2017.10.049)。


**指标和字段怎样理解**

Level3为1,319,138个实例×12,328特征；Level5为473,647条签名×12,328特征。其中205,034条属于小分子签名，其余不能一概算作药物处理。

两个发布的pr_is_lm=1集合相同，都是978个实测基因；Phase I另有gene_info_delta_landmark也是978行，却与这组只重合785个。按ID和字段选择，不能只看行数。

指标解释同Phase II：L3归一化表达、L5响应z，均与CCLE TPM、功能活性分开。


**下载了什么，从哪里下载**

两张大矩阵压缩文件及说明已经在工作区下载校验；本轮只读取了其元数据和细胞覆盖，尚未提取这张大矩阵训练。文件名维度和sig_info行数相互核对，不能把“已下载”写成“已经加入联合模型”。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz) |
| `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz) |
| `GSE92742_Broad_LINCS_sig_info.txt.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz) |
| `GSE92742_Broad_LINCS_cell_info.txt.gz` | [下载文件](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742_Broad_LINCS_cell_info.txt.gz) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

用于增加表达监督或预训练。先审计跨发布的签名/结构重复，保护测试药物及结构组，再将只有表达答案的条件加入部分标签设计。


[回到来源导航](#按来源阅读)


<a id="mixseq"></a>

## 9. MIX-seq：同一扰动下多细胞背景的单细胞RNA

**介绍在哪里，谁产生的数据**

McFarland、Paolella等Broad及合作团队开发的混合细胞单细胞转录方法，研究共同响应和细胞背景特异响应，并分析表达与药敏之间的关系。


[研究介绍与方法](https://doi.org/10.1038/s41467-020-17440-w) · [作者数据发布页](https://doi.org/10.6084/m9.figshare.10298696)


**有哪些细胞系，实际有多少数据**

本地26个表达ZIP共154,710个条形码；按作者normal标记且有DepMap ID的130,648个细胞，实际覆盖169个DepMap ID。完整名单按原ID连接24Q4 Model表取得名称，并保留每个模型实际回收细胞数；未在该版Model表取到名称时继续显示原ID。

这个169是26包中可识别normal细胞的实际联合覆盖，不是所有实验都测了全部169个细胞。具体实验池覆盖见分类文件。

从实际覆盖中举例：A549、PC-3、22Rv1、253J-BV、42-MG-BA、5637、639-V、769-P、786-O、8505C、A-375、A2780。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/mixseq.csv)。


**实验怎样测出来**

把多个细胞系混合培养并施加某一药物/条件，再进行单细胞RNA测序；用SNP遗传信息与参考基因型识别各细胞所属细胞系，得到基因×细胞的稀疏计数和分类/QC。对照也按相应实验池采集。[原始研究](https://doi.org/10.1038/s41467-020-17440-w)和[作者发布](https://doi.org/10.6084/m9.figshare.10298696)。


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| matrix.mtx | 每个基因在每个条形码的RNA计数 |
| classifications.csv | DepMap_ID、细胞分类和cell_quality；normal为作者分类标记 |
| pseudo-bulk | 后续可按细胞系×实验重复汇总的计数，保留重复和对照 |
| all_CL_features.rds中的sens | 作者定义为1−AUC_avg；AUC来自GDSC/PRISM并经其归一化汇总，不是这次RNA测序新测的一条剂量级活性 |

本地包中有Trametinib和DMSO的6h、24h实验，其他包按作者实验说明及补充表定位剂量/时间，不仅凭药名赋值。


**下载了什么，从哪里下载**

26个单药/对照ZIP，每包通常有matrix.mtx、genes.tsv、barcodes.tsv、classifications.csv；另有README、实验/药物补充表和all_CL_features.rds。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `Trametinib_24hr_expt1.zip` | [下载文件](https://ndownloader.figshare.com/files/18716339) |
| `DMSO_24hr_expt1.zip` | [下载文件](https://ndownloader.figshare.com/files/18716354) |
| `README.txt` | [下载文件](https://ndownloader.figshare.com/files/24848579) |
| `Supplementary Tables.xlsx` | [下载文件](https://ndownloader.figshare.com/files/18716510) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

可做跨平台表达复验、不同细胞背景响应和群体异质性研究。先按DepMap ID及实验条件分组，与同批vehicle对照形成可比较表达，再分别追溯其附带功能来源。本轮没有把这些单细胞直接当HMS独立功能样本。


[回到来源导航](#按来源阅读)


<a id="sciplex2"></a>

## 10. sci-Plex2：A549四药剂量实验

**介绍在哪里，谁产生的数据**

Srivatsan等的sci-Plex原始研究中，用A549四药实验验证核标签与多样本表达测量。实验作者产生数据，scPerturb另将其整理成当前H5AD，因此实验出处和下载出处分别标注。


[原始研究](https://doi.org/10.1126/science.aax6234) · [原实验归档](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139944) · [本地下载版本及说明](https://zenodo.org/records/13350497) · [整理脚本](https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/SrivatsanTrapnell2020.py)


**有哪些细胞系，实际有多少数据**

只有**A549**。本地scPerturb 1.4文件有24,262个细胞/核、58,347个特征；药物标签为BMS、Dex、Nutlin、SAHA及control，分别对应BMS-345541、dexamethasone、nutlin-3a、vorinostat及整理后的对照。

从实际覆盖中举例：A549。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/sciplex2.csv)。


**实验怎样测出来**

在独立孔处理A549，作者设计24h、四药、各七个剂量和三重复；用核hash寡核苷酸标记来源，混样后做组合索引单细胞RNA测序，再恢复细胞/核与处理孔的关系。这里的RNA计数是单细胞表达；原文另分析回收细胞数与毒性关系。[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7289078/)。


**指标和字段怎样理解**

| 字段/指标 | 本地实际含义 |
|---|---|
| X | 稀疏RNA计数；一行是细胞/核，一列为特征 |
| cell_line / perturbation / well | 整理后的细胞、处理及孔注释 |
| ncounts / ngenes | 每个细胞的RNA计数与检出特征质量信息，不是活性读数 |
| dose_value | 字符串0、0.1、0.5、1、5、10、50、100；529行缺标签，文件没有dose_unit |
| control | 整理脚本把缺失处理标签设为control；后续须回查原孔/vehicle，不能自动当成已核实的一类DMSO样本 |

**剂量字段已查到来源：**scPerturb代码把top_oligo用下划线分成处理、剂量标签和孔，没有补物理单位。原文Fig.2的Dex轴为0.05—50 μM、Nutlin轴为0.25—250 μM，与本地通用标签序列不同。因此当前保留原始标签，下一步按药物和原实验孔/补充表映射真实剂量，不直接套用sci-Plex3的nM单位。[整理代码](https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/SrivatsanTrapnell2020.py)、[原文Fig.2](https://www.srivatsan-lab.com/pdfs/sci-plex.pdf)。


**下载了什么，从哪里下载**

下载文件为SrivatsanTrapnell2020_sciplex2.h5ad，来自scPerturb Zenodo 1.4；X为表达，obs为细胞和处理注释，var为基因注释。原实验在GSE139944，不把整理版当成另一批独立实验。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `SrivatsanTrapnell2020_sciplex2.h5ad` | [下载文件](https://zenodo.org/api/records/13350497/files/SrivatsanTrapnell2020_sciplex2.h5ad/content) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

适合先验证处理/对照、基因映射与伪bulk流程。剂量和对照注释先按上面的具体来源补足，再做与浓度相关的对接；当前未参与主模型训练。


[回到来源导航](#按来源阅读)


<a id="sciplex3"></a>

## 11. sci-Plex3：三种细胞的大规模药物RNA响应

**介绍在哪里，谁产生的数据**

同一sci-Plex研究的大筛选部分。作者用188种化合物处理A549（肺癌）、K562（髓系白血病）、MCF7（乳腺癌），观察药物模式及单细胞异质性。scPerturb 1.4发布说明标记了对此数据集的修正。


[原始研究](https://doi.org/10.1126/science.aax6234) · [GEO归档](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139944) · [scPerturb 1.4发布](https://zenodo.org/records/13350497)


**有哪些细胞系，实际有多少数据**

实际有注释的三种细胞及行数为：**A549 244,281；K562 173,652；MCF7 344,862**。另36,522行缺细胞/药物注释，共799,317行。

A549包含162,171行24h和82,110行72h；K562、MCF7所列均为24h。72h不是三细胞平衡时间序列。文件有110,983特征，其中人ENSG前缀58,347、鼠ENSMUSG前缀52,636；这是特征列注释，不代表有52,636个鼠细胞。

从实际覆盖中举例：A549、K562、MCF7。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/sciplex3.csv)。


**实验怎样测出来**

不同处理孔的核加样本hash标签，混合后做sci-RNA-seq3，从RNA及hash序列恢复单细胞表达和药物/剂量身份。原研究大筛选设计为四剂量、两重复，并设vehicle；论文分析规模649,340个细胞与本地整理文件799,317行分别使用。[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7289078/)、[整理版本说明](https://zenodo.org/records/13350497)。


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| X / ncounts | RNA计数与每细胞计数总量 |
| cell_line / perturbation | 细胞背景和处理名称；189个非缺失处理名含control |
| dose_value / dose_unit | 本地明确为nM：0、10、100、1000、10000，另有缺失 |
| time / replicate / plate / well | 时间、重复、板、孔，构造对照与伪bulk的依据 |

10、100、1000、10000 nM分别为0.01、0.1、1、10 μM。作者以核回收数量估计的viability和RNA计数是不同量；若用于功能，须说明估计方式及独立计数验证范围，不直接把ncounts称为活性。


**下载了什么，从哪里下载**

下载的是scPerturb 1.4的H5AD，而非原GEO原始文库。处理时先选择人基因、处理重复基因/缺注释，再按实际实验重复汇总。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `SrivatsanTrapnell2020_sciplex3.h5ad` | [下载文件](https://zenodo.org/api/records/13350497/files/SrivatsanTrapnell2020_sciplex3.h5ad/content) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

用于检验跨实验表达预测及群体组成变化。先保留单细胞异质性和同条件对照，分别评价群体响应与细胞状态；功能预测继续以明确的独立测量或经验证估计定义目标。


[回到来源导航](#按来源阅读)


<a id="nci_tpw"></a>

## 12. NCI TPW / GSE116436：多时间点芯片表达

**介绍在哪里，谁产生的数据**

美国NCI/Frederick团队建立的Transcriptional Pharmacodynamics Workbench，用NCI-60模型研究抗癌药物的动态转录反应及其与敏感性的联系。


[研究介绍](https://pubmed.ncbi.nlm.nih.gov/30355619/) · [方法原文](https://pmc.ncbi.nlm.nih.gov/articles/pmid/30355619/) · [GEO发布](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116436) · [TPW工具](https://tpwb.nci.nih.gov/GeneExpressionNCI60/index.html)


**有哪些细胞系，实际有多少数据**

本地GEO样本表头有60个细胞来源标签、7,623个表达样本。包含MCF7、MDA-MB-231、HS-578T、T-47D、MDA-MB-468、A549、HCT-116、PC-3、K-562、U251等；60个完整名称及各自样本数见CSV。

从实际覆盖中举例：MCF7、MDA-MB-231、A549、HCT-116、PC-3、K-562、HT29、786-0、A498、ACHN、BT-549、CAKI-1。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/nci_tpw.csv)。


**实验怎样测出来**

细胞给药后在2、6、24h提取RNA，标记并与Affymetrix U133A平台探针杂交，洗涤、扫描得到探针信号。原方法以RMA处理背景和归一化；TPW工具再以处理与相应零浓度的log2表达差描述响应。[原研究方法](https://pmc.ncbi.nlm.nih.gov/articles/pmid/30355619/)。


**指标和字段怎样理解**

| 指标/字段 | 意义 |
|---|---|
| GEO处理后探针表达 | 芯片探针×样本数值，需按平台注释对应基因；不是RNA测序UMI |
| log2表达差 | 同细胞/时间下处理值减相应对照值，作为后续可构造的响应 |
| 样本标题/来源 | 如786-0_5-Azacytidine_1000nM_24h，解析细胞、药物、nM浓度和时间 |
| 0nM | 相关零剂量参照，仍需核对细胞、时间与实验批次 |

2h/6h/24h分别2,540/2,551/2,532列；零剂量2,554列。平台标识GPL571。原研究设计含15药和两种非零浓度；本地样本依标题逐条定位。


**下载了什么，从哪里下载**

下载GSE116436_series_matrix.txt.gz，是GEO处理后表达及样本说明文件。细胞名单来自该文件source_name_ch1表头，而非另找一张通用NCI-60名单替代。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `GSE116436_series_matrix.txt.gz` | [下载文件](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116436/matrix/GSE116436_series_matrix.txt.gz) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

可验证时间变化和外部表达泛化。下一步按GPL571映射探针，再按细胞、药物、剂量、时间及对照构造表达响应；当前已读表头与条件，尚未完成全矩阵的基因级对接。


[回到来源导航](#按来源阅读)


<a id="drugseq"></a>

## 13. Novartis DRUG-seq U2OS MoABox：每孔转录组

**介绍在哪里，谁产生的数据**

Novartis公开高通量DRUG-seq MoABox表达资料，帮助研究化合物作用机制和开发分析方法。作者介绍4,343种小分子；本地比较成员表有4,344个cmpd_sample_id，样品ID/对照和唯一药物分别统计。


[作者项目、测量及对象说明](https://github.com/Novartis/DRUG-seq/blob/main/data/Novartis_drugseq_U2OS_MoABox/README.md) · [固定数据发布](https://zenodo.org/records/14291446)


**有哪些细胞系，实际有多少数据**

仅一种细胞：**U2OS / U-2 OS，骨肉瘤来源**。52批、每批3块384孔板，共156块板、59,904个孔列。作者标注3,276个EMPTY孔，56,628个样本参与后续处理。59904是孔列数，不是59904种细胞或单细胞数。

从实际覆盖中举例：U2OS。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/drugseq.csv)。


**实验怎样测出来**

在U2OS培养孔中进行药物处理，24h后用带样本条形码和UMI的高通量RNA测量流程得到每孔基因计数。重复板保留三个实验重复；它是孔内群体RNA测量，sample ID对应孔。作者README还描述DMSO参照孔、阳性对照及EMPTY孔。[实验说明](https://github.com/Novartis/DRUG-seq/blob/main/data/Novartis_drugseq_U2OS_MoABox/README.md)。


**指标和字段怎样理解**

| 对象/字段 | 意义 |
|---|---|
| Exp / UMI.counts | 每板59,594特征×384孔的RNA分子计数 |
| Annotation | 对应孔、样品及处理注释 |
| external_biosample_id / plate_well | 具体测量样本、板与孔 |
| cmpd_sample_id / concentration / unit | 药物样品、浓度、单位 |
| comparisons_info_public | 比较设计，处理与参照的关系；行数不是新的实验样本数 |

作者常规剂量0.01、0.1、1、10 μM，时间24h；DMSO为参考。下载的两个RData没有单独的细胞活性实测标签；作者另外发布的差异表达/PCA等与UMI原计数分别处理。


**下载了什么，从哪里下载**

Exp_gzip.RData保存Exp及板级计数/注释；comparisons_metadata.RData保存比较设计及56,628条样本成员记录。本地已实际读对象与维度。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `Exp_gzip.RData` | [下载文件](https://zenodo.org/api/records/14291446/files/Exp_gzip.RData/content) |
| `comparisons_metadata.RData` | [下载文件](https://zenodo.org/api/records/14291446/files/comparisons_metadata.RData/content) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

适合补充大量药物的表达监督及结构泛化测试。按cmpd_sample_id连接化合物结构，按比较成员表选同板/批次参照，形成24h表达目标。只覆盖U2OS这一背景，后续迁移评价与本轮五细胞分开。


[回到来源导航](#按来源阅读)


<a id="tahoe_shard"></a>

## 14. Tahoe-100M：大规模单细胞药物响应，当前一个分片

**介绍在哪里，谁产生的数据**

Tahoe官方在Hugging Face公开、由Mosaic高通量平台生成的药物单细胞图谱，目标是描述不同细胞背景中的药物分子响应。


[研究介绍、字段和下载结构](https://huggingface.co/datasets/tahoebio/Tahoe-100M) · [官方README](https://huggingface.co/datasets/tahoebio/Tahoe-100M/blob/main/README.md) · [文件目录](https://huggingface.co/datasets/tahoebio/Tahoe-100M/tree/main)


**有哪些细胞系，实际有多少数据**

作者图谱介绍约50种癌细胞背景和约100M转录组。本地只下载第1个表达分片：**28,225个细胞，50个实际cell_line_id，66个sample、66个处理名称**。50个完整名称是从表达分片的Cellosaurus ID连接细胞注释表得到的。

细胞注释表本身1,000行，去重后有102个身份，许多行描述不同driver注释；所以既不能把1000行称为1000个细胞系，也不能默认102个身份都出现在当前表达分片。

从实际覆盖中举例：A549、A-172、A-427、A498、AN3 CA、AsPC-1、BT-474、C-33 A、C32、CFPAC-1、CHP-212、COLO 205。

[下载全部名单、原ID和每个细胞的记录数（CSV）](../data/catalog/tahoe_shard.csv)。


**实验怎样测出来**

在多细胞背景的高通量培养/干预设计中给药，再以单细胞RNA测序和身份注释取得各细胞的稀疏基因计数。官方数据卡说明mixed-cell spheroid、96孔处理板及plate标识。[官方数据卡及字段](https://huggingface.co/datasets/tahoebio/Tahoe-100M/blob/main/README.md)。本地sample元数据含浓度组合，但没有独立逐样本处理时长列，本次不从其他实验补写统一时间。


**指标和字段怎样理解**

| 字段/指标 | 意义 |
|---|---|
| genes / expressions | 对齐的基因token与原RNA计数列表；首项CLS/marker需按说明处理 |
| cell_line_id | 实际细胞身份，使用Cellosaurus ID；表中另有DepMap映射 |
| sample / plate | 生物样本/处理重复及板来源 |
| drug / canonical_smiles | 处理与结构；DMSO_TF为vehicle标识，结合plate找参照 |
| drugname_drugconc | 样本说明中的药物—浓度—单位组合 |
| targets / moa | 药物注释，部分由模型辅助整理，不是独立RNA或功能实验测量 |

当前主要数值为单细胞RNA counts；“每细胞RNA数”与“细胞活性”不同。379条药物说明、1,344条sample说明的覆盖也不等于当前66个样本/处理名称。


**下载了什么，从哪里下载**

四个Parquet文件：train-00000-of-03388为一个表达分片，另三个为药物、样本、细胞注释。全图谱、元数据表和本地一个表达分片是三个不同范围。


| 本地文件名 | 本次下载记录的文件直链 |
|---|---|
| `train-00000-of-03388.parquet` | [下载文件](https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/data/train-00000-of-03388.parquet) |
| `drug_metadata.parquet` | [下载文件](https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/drug_metadata.parquet) |
| `sample_metadata.parquet` | [下载文件](https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/sample_metadata.parquet) |
| `cell_line_metadata.parquet` | [下载文件](https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/cell_line_metadata.parquet) |

以上是该来源的主要文件。其余说明、原孔或支持文件的准确URL、版本和用途见[完整下载清单](../data/download_manifest.csv)。原始大文件在工作区3_MVP/data/raw，已解压表在3_MVP/data/unpacked；仓库保存当前训练所需文件及这份来源/覆盖记录。


**为什么使用，怎样接入，实际做到哪里**

先验证格式、基因token及实际条件；再按需要选择分片和对照用于表达预训练/状态分析。当前完成本地分片及实际身份核查，没有将图谱宣传规模计入本轮训练。


[回到来源导航](#按来源阅读)


## 本次核查和后续使用的共同规则

- 细胞名单从实际RNA矩阵行、实际签名、实际药物孔/浓度表、实际单细胞注释产生；说明表更广的覆盖另外报告。提取代码是[src/catalog_data.py](../src/catalog_data.py)，每项计数依据见[summary.json](../data/catalog/summary.json)。
- 测量值必须保留定义：TPM、z响应、UMI、原始强度、相对计数、条形码变化、曲线AUC各自代表不同量。身份、结构、实际剂量和时间查清楚后，再决定同条件、部分可比或只有单类标签的学习方式。
- 本次已核对可读的Sanger方法、PRISM项目、Novartis说明、scPerturb发布/代码、Tahoe数据卡和原论文。部分DepMap/GEO入口返回浏览器验证，HMS/CTRP入口在本次工具访问中未正常返回正文；链接来源仍可追到固定下载清单及当时实际保存文件，测量说明同时核对作者论文、README或Excel Meta。没有把网页暂时不可读解释成数据不存在。
- 新补充的覆盖表只说明数据资源及条件，不增加本轮训练样本；模型、训练划分、权重和指标没有因本次文档整理改变。

[首页](../README.md) · ← [研究框架](01_framework.md) · [下一步看：具体怎样对接](03_integration.md) →
