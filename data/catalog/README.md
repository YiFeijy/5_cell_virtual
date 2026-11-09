# 实际细胞覆盖清单

由固定本地文件提取，2026-11-09。主说明见[数据页](../../docs/02_data.md)。每个CSV是一个实际文件范围，命名/ID原样保留，行数不是所有来源之间去重后的细胞系数量。

| 清单 | 范围及主要计数 |
|---|---|
| [ccle_24q4.csv](ccle_24q4.csv) | 1673条身份/矩阵行记录；原单位和范围见summary.json |
| [lincs_gse70138.csv](lincs_gse70138.csv) | 41条身份/矩阵行记录；原单位和范围见summary.json |
| [lincs_gse92742.csv](lincs_gse92742.csv) | 76条身份/矩阵行记录；原单位和范围见summary.json |
| [hms.csv](hms.csv) | 6条身份/矩阵行记录；原单位和范围见summary.json |
| [ctrpv2_2015.csv](ctrpv2_2015.csv) | 887条身份/矩阵行记录；原单位和范围见summary.json |
| [ctrpv2_2016.csv](ctrpv2_2016.csv) | 835条身份/矩阵行记录；原单位和范围见summary.json |
| [gdsc1.csv](gdsc1.csv) | 970条身份/矩阵行记录；原单位和范围见summary.json |
| [gdsc2.csv](gdsc2.csv) | 969条身份/矩阵行记录；原单位和范围见summary.json |
| [prism_primary.csv](prism_primary.csv) | 578条身份/矩阵行记录；原单位和范围见summary.json |
| [prism_secondary.csv](prism_secondary.csv) | 489条身份/矩阵行记录；原单位和范围见summary.json |
| [mixseq.csv](mixseq.csv) | 169条身份/矩阵行记录；原单位和范围见summary.json |
| [sciplex2.csv](sciplex2.csv) | 1条身份/矩阵行记录；原单位和范围见summary.json |
| [sciplex3.csv](sciplex3.csv) | 3条身份/矩阵行记录；原单位和范围见summary.json |
| [nci_tpw.csv](nci_tpw.csv) | 60条身份/矩阵行记录；原单位和范围见summary.json |
| [drugseq.csv](drugseq.csv) | 1条身份/矩阵行记录；原单位和范围见summary.json |
| [tahoe_shard.csv](tahoe_shard.csv) | 50条身份/矩阵行记录；原单位和范围见summary.json |

## 怎样读

- `cell_name`是源名称，或依据精确ID连接到固定Model说明的名称；MIX-seq取名来源另有字段。
- `model_id`、`sanger_model_id`、`cell_id`、`Cell_ID_Cellosaur`保留相应来源的身份体系。
- `raw_well_rows`、`dose_rows`、`all_signatures`、`cell_or_nucleus_rows`、`expression_samples`分别是孔、浓度记录、汇总签名、单细胞/核及表达样本数，不混成一种样本数。
- PRISM保存实际矩阵行及STR问题；失败行不自动变成一个认证细胞身份。
- sci-Plex3的36,522条缺注释行在summary中单列，不凭表达相似性补细胞名称。
- Tahoe名单仅列已下载表达分片实际出现的50个身份，注释表的102个不同身份另记。

复现命令（需要完整原始文件，`--source-root`指包含data/raw及data/unpacked的工作区）：

```bash
python src/catalog_data.py --source-root ..
```

数字摘要：[summary.json](summary.json)。介绍页与文件出处索引：[source_pages.csv](source_pages.csv)。完整下载文件URL仍使用原下载清单，不再复制大文件。

## PharmacoSet整理对象的补充核查

[pharmacoset.json](pharmacoset.json)记录直接readRDS的对象维度、与2015实际细胞ID集合的核对，以及两条tipifarnib源ID合并的情况。它复用同一887细胞集合，不另加一张重复名单。R中可以直接检查：

```r
x <- readRDS("data/raw/CTRPv2/CTRPv2_PharmacoSet.rds")
dim(attr(x, "cell"))
dim(attr(x, "drug"))
dim(attr(x, "sensitivity")$raw)
```

下载该RDS后执行；其他原始覆盖表由Python目录脚本生成。LINCS清单中small_molecule_signatures大于0的行才是小分子子集，all_signatures包括其他扰动和对照。
