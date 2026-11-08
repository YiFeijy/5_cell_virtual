# 07 检查、训练与数据重建

先理解[研究步骤](../README.md)，本页说明如何检查同一条证据链。建议Python 3.12，本轮实际包版本保存在[环境记录](../results/environment.json)。

## 1. 不训练，先核对已有数据和模型

```bash
git clone https://github.com/YiFeijy/5_cell_virtual.git
cd 5_cell_virtual
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python src/check.py
python src/inspect_one.py
```

`check.py`检查文件校验、样本/数组顺序、跨库板权重、划分隔离及已发布功能预测对应的指标。`inspect_one.py`读取一个已有权重，重算示例，不训练。示例属于未见药物第1折测试集，种子20261107；其他折的样本需要先生成相应权重。

## 2. 输入与答案怎样读取

```python
import numpy as np
import pandas as pd

samples = pd.read_csv('data/training/samples.csv')
with np.load('data/training/arrays.npz', allow_pickle=False) as a:
    i = 311
    background = a['baseline'][a['cell_index'][i]]
    fingerprint = a['drug_morgan'][a['drug_index'][i]]
    log_dose = a['log10_dose_uM'][i]
    expression_answer = a['expression_z24'][i]
    count_answer = a['count72'][i]
    print(samples.loc[i, ['cell_id', 'pert_iname', 'dose_uM']])
    print(background.shape, fingerprint.shape, log_dose,
          expression_answer.shape, count_answer)
```

背景、结构、剂量是输入；两个实测答案用于监督和评价。缩放只拟合训练部分，各次运行的缩放参数保存在其`preprocessing.npz`中。

## 3. 用已发布训练数据重跑模型

在独立克隆中运行，以便保留已有结果供比较：

```bash
python src/train.py --workers 6
python src/baselines.py
python src/analyze.py
python src/information_diagnostic.py
python src/plots.py
```

顺序对应：144次神经比较→简单/岭回归参照→测试与排序→额外实测表达诊断→图。训练跳过仓库附带的那次完整运行，其余运行重新生成。144次历史记录已合并到`results/training_history.csv.gz`，不会被误判成144个可跳过的完整模型。

冻结配置核对数组、源码及torch版本。仓库目录整理只改变训练代码的根目录/输入路径，网络、损失、优化与划分计算保持一致；原训练源码哈希与迁移验证保存在[provenance.json](../results/provenance.json)。本次整理没有新训练。新方法或数据修改应另建实验版本，保留原结果并冻结新的评价设计。

## 4. 从原始文件核算数据

本轮数据链为：**HMS原孔核算→已核对的条件/身份表→Level5表达提取→Level3实例与对照QC→训练数组**。

HMS三批Excel及LINCS说明表随仓库保留；大矩阵从清单获取。源文件下载支持续传、校验和只列清单：

```bash
python -m pip install -r requirements-data.txt
python src/download_data.py --dataset CCLE_DepMap_24Q4 --dataset LINCS_GSE70138 --dataset HMS_LINCS --dry-run
python src/download_data.py --dataset CCLE_DepMap_24Q4 --dataset LINCS_GSE70138 --dataset HMS_LINCS --workers 4
python src/unpack_lincs.py
python src/aggregate_hms.py
python src/extract_reference.py
python src/prepare.py
```

三个来源原下载文件约18.1 GiB，解压另占空间。下载清单中还有其他来源；`--dataset`限定本次范围。

`aggregate_hms.py`重算每批技术孔和三批均值，并与保留的功能标签逐条核对。`extract_reference.py`按已核查的签名ID、基因ID从原始Level5重新提取当前3,410条表达候选；`prepare.py`再读取原始Level3实例和对照，重建3,364条训练观测及划分。

`data/reference/`中的身份/结构/条件对应表是前期核查产物，随仓库保留作为本轮固定输入定义；它们与原始文件、训练集分别标识。MCF10A和受限结构药物的处理原因见[数据说明](02_data.md)。这条重建链围绕当前冻结条件，不声称恢复了所有公开来源的完整整合。

第三方数据和文献保留来源，使用遵循原发布方条款。
