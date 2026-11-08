# 药物分子响应与功能数据调研

本轮完成公开数据与相关文献的广泛检索，并对Broad–HMS、LINCS元数据、MIX-seq原始表达和部分PRISM、Tahoe、DepMap元数据进行实际读取。原有《初步方案.md》和1.png保持原样。

建议按以下顺序阅读：

1. [数据可行性审查](data_audit.md)：读数定义、优先数据、实际下载核验和当前缺口。
2. [实际交集与匹配规模](overlap_report.md)：区分候选条件数、完整输入元数据和矩阵可用性。
3. [已有研究怎样对接](literature_review.md)：18项核心文献与关键作者代码。
4. [十二类对接方法](integration_methods.md)：实验配对、近似浓度、曲线插值、非配对联合学习及时间关联。
5. [数据驱动的验证顺序](experiment_plan.md)：下一阶段的任务与完成条件。

机器可读资源：[50项资源清单](dataset_catalog.csv)、[113条细胞和药物映射](data_mapping.csv)、[下载来源](download_sources.csv)、[文件大小与哈希](audit/downloaded_file_manifest.csv)。资源家族、整合入口和不同版本有重叠，不能把条目数当作独立配对数据集数。

已读取的数据与作者文件位于`metadata/`及`evidence/`；统计、候选配对和MIX-seq试读伪bulk位于`audit/`。没有训练预测模型，也没有把跨实验记录当作同一个细胞的真实连续轨迹。

重算顺序：

```bash
python3 数据可行性调研_20261105/scripts/fetch_metadata.py --include-expression
python3 数据可行性调研_20261105/scripts/build_catalog.py
python3 数据可行性调研_20261105/scripts/audit_metadata.py
Rscript 数据可行性调研_20261105/scripts/audit_mixseq.R
python3 数据可行性调研_20261105/scripts/audit_mixseq.py
python3 数据可行性调研_20261105/scripts/summarize_audit.py
```

首次重新下载时需能访问相应公共网站；来源返回HTML时，不得据HTTP成功状态判断数据下载成功。`fetch_metadata.py`对数据文件检查基本格式并记录结果；哈希用于固定本次本地内容，不等同于所有来源均提供了官方SHA256。
