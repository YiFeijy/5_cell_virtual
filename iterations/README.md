# 从零手动复现路线图与迭代台账（分支 `step-by-step-repro`）

> **本分支目标**：面向“**自己亲手一步一步从头复现**”，把整个抗癌药物筛选虚拟细胞项目（工作一）的演化过程拆成**颗粒度清晰、每次迭代独立提交**的可运行步骤。
>
> **铁律约定（每一次迭代单独 1 次 Git Commit）**：
> 每次提交必须把以下 6 件事写透并落盘可核验产物：
> 1. **思路（Why）**：这一步要回答什么科学/工程问题？
> 2. **数据（Data）**：输入来自哪个原始文件、字段/单位是什么、做了哪些筛选/对齐/防泄漏？
> 3. **模型与计算策略（Method）**：具体网络结构、损失函数、划分协议与对照基线是什么？
> 4. **代码（Code）**：一份清晰可读、带中文注释、支持本地直接手动运行的独立脚本。
> 5. **改动（Diff）**：相比上一步改了什么？为什么要这么改？
> 6. **实测结果与诊断（Results → Next）**：真实运行后的完整数值（含正结果与负结果），以及由结果直接引出的下一步迭代方向。

---

## 一、总体研究思路（为什么要“先做 v0，再逐代演进”）

我们想回答的终极问题是：
> **给癌细胞施加某种药物和剂量后，既会改变细胞的基因表达状态，也会改变细胞的存活/增长数量。能否先从最干净的配对数据出发，检验“表达 + 细胞数联合建模”是否真的有助于预测未见新药与筛选候选药物？**

但这个领域最容易踩的两个坑是：
1. **一上来就堆复杂模型和大杂烩数据**，出了分也不知道是哪部分数据的功劳、有没有化学泄漏；
2. **跳过简单基线**（如 `cell×dose` 均值和线性岭回归），误把看似不错的数字当成深度学习优势。

因此，本分支严格按照**“先小而干净的最初闭环 → 暴露真问题 → 按证据逐次迭代升级”**的科研节奏推进：

```text
iter-00  盘点已下载原始数据：每个文件长什么样、字段单位、数据陷阱（MCF10A 缺 RNA、-666 哨兵等）
   │
   ▼
iter-01  最初版本 (v0) 数据预处理：只用来源最干净的同项目平行数据（CCLE + LINCS GSE70138 + HMS）
         36,912 原孔 → 4,140 库板条件 → 3,364 条质控观测（5 细胞 × 106 药 × 24h 表达 × 72h 核计数）
   │
   ▼
iter-02  最初版本 (v0) 第一部分计算策略与模型：
         简单参照（均值 / 岭回归）vs 单任务 vs 联合（共享 / 私有 / 瓶颈）+ 事后实测表达信息诊断
         ──► 真实暴露 v0 的瓶颈：106 药太少，新药连表达都打不过均值；但实测 24h 表达可极准预测 72h（MAE 0.125）！
   │
   ▼  （后续迭代将在本分支依次单独提交）
iter-03  v1 迭代：接入 LINCS 2020（~3 万药）+ CTRPv2（多谱系 496 药）解决小样本化学泛化瓶颈
iter-04  v2 迭代：反泄漏审计（互变异构体 / 手性）与严格单变量因果归因（no_lincs / shuffled_lincs）
iter-05  v3 迭代：双轨动态状态空间模型（dyn_state_only）+ 19k 受控 CCLE 药理面板 + 生长率校准与 3h/24h 因果同化
```

---

## 二、本分支迭代提交总账

| 迭代目录 | 对应步骤 | 核心内容 | 手动运行入口 | 当前状态 |
| :--- | :--- | :--- | :--- | :--- |
| [00_raw_data_audit](00_raw_data_audit/README.md) | `iter-00` | 已下载原始数据全貌、v0 三大核心库深度剖析、扩展库规划 | `python iterations/00_raw_data_audit/inspect_raw_data.py` | 已完成并提交 |
| [01_v0_preprocessing](01_v0_preprocessing/README.md) | `iter-01` | 最初版本 (v0) 数据清洗、HMS 原孔重算、L3/L5 对齐 QC、特征与 5 折划分 | `python iterations/01_v0_preprocessing/run_preprocess_v0.py` | 已完成并提交 |
| [02_v0_baseline_and_joint_model](02_v0_baseline_and_joint_model/README.md) | `iter-02` | 最初版本 (v0) 计算策略、基线与 8 组神经网络代码、完整复现结果与问题诊断 | `python iterations/02_v0_baseline_and_joint_model/run_v0_experiment.py --full` | 已完成并提交 |
| `03_v1_...`（待后续提交） | `iter-03+` | 扩大药物规模（LINCS 2020 + CTRPv2）→ 反泄漏 → 动态状态建模 | 待后续按轮次独立提交 | 规划中 |

---

## 三、本地手动复现环境准备

1. **克隆本分支**：
   ```bash
   git clone -b step-by-step-repro https://github.com/YiFeijy/5_cell_virtual.git
   cd 5_cell_virtual
   ```
2. **安装依赖**（建议 Python 3.10–3.12）：
   ```bash
   python -m pip install -r requirements.txt
   # 若需直接读取 HMS 原始 .xls Excel，请确保安装 xlrd>=2.0.1：
   python -m pip install xlrd
   ```
3. **按顺序亲手运行前三步**：
   ```bash
   # 第 0 步：先看清楚下载的原始数据长什么样
   python iterations/00_raw_data_audit/inspect_raw_data.py

   # 第 1 步：亲手复现最初版本 (v0) 的完整数据预处理并核对数值一致
   python iterations/01_v0_preprocessing/run_preprocess_v0.py

   # 第 2 步：亲手跑一遍最初版本 (v0) 的基线、神经网络、排序与信息诊断
   #   快速试跑（仅第 0 折 + 1 个种子，约 10 秒）：
   python iterations/02_v0_baseline_and_joint_model/run_v0_experiment.py --quick
   #   完整复现（5 折新药 + 高剂量外推 + 3 种子 + 基线 + Bootstrap CI + 信息诊断，约 1-2 分钟）：
   python iterations/02_v0_baseline_and_joint_model/run_v0_experiment.py --full --workers 6
   ```

> **说明**：本分支不修改 `main` 上原有已冻结的 `src/`、`data/training/` 和 `results/` 文件；所有分步教学与迭代脚本及复现输出均集中在 `iterations/` 目录下，方便逐 commit 对比阅读和独立运行。
