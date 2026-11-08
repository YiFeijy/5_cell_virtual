# 原始数据下载与用途

最新模型进展见[2026-11-07初步训练与探索](../第一部分_初步实验_20261107/README.md)：已完成144个神经训练任务及外层探索评价，并保存基线/排序/实测表达诊断。下文“尚未训练”及v3准备包中的零训练记录描述更早的数据准备阶段。

整体研究设计见[项目总设计](../项目总设计_20261106/README.md)。最新实验准备产物见[HMS v3说明](prepared_20261106/HMS_v3_实验准备/README.md)，包括实际数组、来源、权重、结构与重复审计及划分草案；生物学入组与最终训练发布尚未冻结。

当前对接进度以[数据内容 处理步骤与设计原因](investigation_20261106/当前对接说明_数据内容处理步骤与设计原因.md)为入口。HMS已有[同库板来源核对的v2产物](integration_20261106/v2_库板来源核对/README.md)；[MCF10A缺行的原因与补充来源](investigation_20261106/MCF10A为何缺行_原因追查.md)已另行追查。

先了解数据本身，可以阅读[数据集说明 来源与测量内容](数据集说明_来源与测量内容.md)。它逐项介绍原研究目的、产生数据的机构、测量方法、观测单位、字段与下载来源。

本批下载队列的78个文件已全部完成校验。对接设计、操作顺序及本轮试对接结果见[对接方案与操作步骤](integration_20261106/对接方案与操作步骤.md)。

本目录保存数据收集、下载校验及实际内容检查结果。11月6日已读取部分实际矩阵，核查同一LJP项目内表达与功能读数的条件对应；尚未构建跨库训练集，也未训练模型。

历史v1结果见[数据检查与可行性分析](analysis_20261106/数据检查与可行性分析.md)。其中[3,168个输入与三个终点可用的条件](analysis_20261106/HMS_3168条件_输入与三个终点可用.csv)是早期可用性清单；后续请以v2来源核对、v3准备包和各自QC状态为准。

| 数据目录 | 下载的是什么 | 数据用来做什么 |
|---|---|---|
| `raw/CCLE_DepMap_24Q4/` | 未给药细胞系的基线RNA表达矩阵，以及细胞和RNA记录说明 | 描述给药前的细胞背景 |
| `raw/LINCS_GSE92742/` | LINCS Phase I的Level3表达、Level5变化签名及完整实验说明 | Level3查看处理与对照的分子状态；Level5描述药物引起的表达变化 |
| `raw/LINCS_GSE70138/` | LINCS Phase II的对应表达、变化签名及实验说明 | 增加药后表达数据和实验条件覆盖；保留不同数据版本 |
| `raw/CTRPv2/` | 原始发布ZIP，以及单独保存的PharmacoSet整理版 | 获得药物浓度级相对活性、药敏曲线和相关实验信息 |
| `raw/GDSC/` | GDSC1和GDSC2原始浓度实验包，以及拟合结果Excel | 原始包提供功能assay读数；拟合表提供IC50、AUC等药敏汇总 |
| `raw/PRISM_19Q4/` | Primary及Secondary相对细胞丰度矩阵、处理说明、细胞说明和曲线参数 | 查看给药后的细胞丰度变化与药物敏感性；其中Secondary包含浓度系列 |
| `raw/HMS_LINCS/` | 平行实验的处理前后细胞计数、相对细胞数和相对生长率 | 直接观察细胞积累或生长受到的影响 |
| `raw/MIXseq/` | 作者发布的单药及对照单细胞表达包、实验说明和敏感性资料 | 观察不同细胞系及同一细胞系内部的响应差异 |
| `raw/sciPlex/` | sci-Plex 2和3的完整标准化单细胞表达文件 | 观察药物、浓度与单细胞分子响应的关系 |
| `raw/NCI_TPW/` | NCI药物处理后多个时间点的表达矩阵 | 观察分子响应随时间如何变化 |
| `raw/Novartis_DRUGseq/` | U2OS药物处理后的表达RData对象及实验条件对象 | 补充全转录组药物响应数据 |
| `raw/Tahoe_100M/` | 一个完整官方表达分片，以及细胞、药物、样本说明表 | 查看Tahoe单细胞数据内容；本批不是Tahoe全库 |

**表达数据描述分子变化。功能assay描述特定测量下的活性、生长或敏感性。文件下载完整之后，再依据实验说明确定它们能否关联。**

药物结构随各数据集说明表一起保存：LINCS的`pert_info`、PRISM的`treatment-info`和Tahoe的`drug_metadata`。这些文件说明药物身份与分子结构；剂量、时间及重复信息在相应实验说明表中。

## 当前进度和文件位置

- [逐文件用途与实时下载状态](数据用途与下载状态.md)
- [完整下载队列与来源地址](download_manifest.csv)
- [机器可读状态](download_status.json)
- [已下载文件的内容检查清单](下载文件内容清单.csv)
- [解压后的CTRP和GDSC文件清单](unpacked/manifest.json)：原始ZIP也继续保留。
- [已解压LINCS矩阵的位置和维度](unpacked/LINCS_matrix_contents.json)：已完成矩阵可直接读取，不需要重新下载。
- `descriptions/`：从原包取出的README、字段字典和GDSC原始读数说明。
- [下载总日志](download_logs/manager.log)
- [11月6日续传日志](download_logs/manager_20261106.log)

`.part`表示尚未完成，不能作为完整数据使用。只有下载结束且通过校验的文件才标记为“已完成”。来源提供校验码时核验官方MD5或SHA512；ZIP还检查CRC。全部保留来源原文件，不在下载过程中做数据合并。

下载由后台tmux会话`mvp_data_download_20261106`继续执行。每个文件支持断点续传，连接中断后会重新读取已下载长度并重试；状态每5秒更新。已完成且未变化的文件复用已有校验记录。重新启动同一个下载程序会继续未完成文件：

```bash
python3 data/scripts/download_all.py --workers 6
```

只允许一个下载管理程序运行，避免同时写入同一文件。每次下载的地址、文件用途、范围和校验结果均保存在清单中。

## 来源

[DepMap固定发布](https://api.figshare.com/v2/articles/27993248)、[LINCS Phase I](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)、[LINCS Phase II](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)、[CTRP维护者下载配置](https://github.com/BHKLAB-DataProcessing/CTRPv2-Pharmacoset_Snakemake/blob/main/config/pipeline.yaml)、[GDSC官方下载页](https://cellmodelpassports.sanger.ac.uk/downloads/gdsc)、[PRISM官方发布](https://depmap.org/repurposing/)、[HMS数据目录](https://lincs.hms.harvard.edu/db/datasets/)、[MIX-seq原始数据](https://doi.org/10.6084/m9.figshare.10298696)、[sci-Plex标准化发布](https://zenodo.org/records/13350497)、[NCI时间序列GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116436)、[Novartis官方归档](https://zenodo.org/records/14291446)、[Tahoe官方发布](https://huggingface.co/datasets/tahoebio/Tahoe-100M)。
