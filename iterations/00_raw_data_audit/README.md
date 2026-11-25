# `iter-00`：总体思路与已下载原始数据全貌剖析

> **本提交目标**：在写任何预处理或模型代码之前，先把本地下载好的全部原始数据翻开看清楚——**有什么文件、多大、每行代表什么、字段和单位是什么、有哪些坑**，再明确最初版本（`v0`）为什么只从其中最干净的 3 个数据集起步。

---

## 1. 本步思路（为什么必须先体检原始数据）

- 公开生物医学数据库存在大量容易踩坑的隐性差异：**细胞系命名不统一、同一批实验分多个库板、不同平台测量的物理量完全不同（核计数 vs ATP 发光 vs 混合条形码丰度）、表达矩阵里混有缺失哨兵值 `-666`、部分药物 SMILES 标为 `restricted`、某些细胞在模型表里有身份但实际没有 RNA 行**。
- 如果不先写脚本把原始文件的真实维度、数值范围和缺行原因查清楚，后面预处理写出来的数字就成了“黑箱”。
- 本步只做**只读体检与盘点**，不建模型、不改原始数据。

---

## 2. 已下载数据总览（本地共 17 个目录 / 100 个文件 / 132.76 GB）

以下统计由 [inspect_raw_data.py](inspect_raw_data.py) 实扫 `/home/data2/jiangyun/5_cell_virtual/3_MVP/data/raw` 自动生成，完整清单见 [raw_data_inventory.csv](raw_data_inventory.csv) 与 [raw_data_profile.json](raw_data_profile.json)：

### 2.1 最初版本（`v0`）直接使用的 3 个核心数据集

| 数据集 | 文件数 / 大小 | 格式 | 在 `v0` 中的角色 |
| :--- | :---: | :--- | :--- |
| **`CCLE_DepMap_24Q4`** | 4 文件 / 484.10 MB | `.csv` | **输入（细胞背景）**：给药前细胞的基线 RNA-seq 表达 `log2(TPM+1)` |
| **`LINCS_GSE70138`**（Phase II） | 8 文件 / 17.62 GB（解压 GCTX 另占 21.32 GB） | `.txt.gz`, `.gctx.gz` | **输入（药物结构/条件）+ 标签（24h 表达）**：978 个实测 landmark 基因的相对响应 z-score |
| **`HMS_LINCS`**（Joint Project） | 3 文件 / 14.09 MB | `.xls` | **标签（72h 功能）**：同项目平行实验的 72h 成像核计数（相对细胞数 `Nt/Nc` 与相对净增长） |

### 2.2 本地已下载、留待后续迭代（`v1 → v3` 及外部盲测）的数据集

| 数据集 | 大小 | 后续迭代用途 |
| :--- | :---: | :--- |
| **`LINCS_2020`** | 33.52 GB | `v1–v3` 核心化学泛化杠杆：`720,216` 条小分子 `trt_cp` 签名（`510,270` 通过 `qc_pass`）、`230` 细胞、`34,419` 药物，含 24h(55.9万) / 6h(13.8万) / 3h(1.2万) / 48h(4,372) / 72h(796) |
| **`LINCS_GSE92742`**（Phase I） | 68.71 GB | L1000 Phase I 原始归档 |
| **`CTRPv2`** | 581.24 MB | `v1–v3` 多谱系 72h 实测 ATP 药敏训练与筛选排序评测（`545` 处理条目中含 `49` 组合、留 `496` 单药，`1,107` 细胞元数据） |
| **`CellModelPassports`** | 182.58 KB | `v3` 训练期细胞固有 72h 对照倍增倍数辅助校准（`1,041` 行 / `964` 模型，含 `seeding_density` 与 `day4_day1_ratio`） |
| **`GDSC`（1 & 2）** | 230.55 MB | 外部 72h 跨测定药敏验证（GDSC1 染料 / GDSC2 ATP） |
| **`PRISM_19Q4` & `PRISM_24Q2`** | 383.65 MB / 1.38 GB | 外部 120h 混合培养条形码药敏盲测与跨时长标度核验 |
| **`NCI_TPW`** | 855.73 MB | 外部跨平台（微阵列芯片）60 细胞 × 15 药 × 2/6/24h 时间演化盲测 |
| **`Novartis_DRUGseq`** | 2.54 GB | 外部 U2OS × 4,343 药跨平台 MoA 检索验证（不含独立活性标签） |
| **`MIXseq` / `sciPlex` / `Tahoe_100M` / `DTP_SingleCell`** | 1.75 / 2.49 / 2.20 / 0.07 GB | 单细胞异质性、时序验证及后续工作二（群体组成解耦与耐药） |
| `HMS_Hafner71_GR` | 0 B | 空占位目录 |

---

## 3. 核心 `v0` 三大原始数据到底长什么样子？（实测剖析）

### 3.1 `HMS_LINCS`（`Screen20237 / 20238 / 20239` 原孔 Excel）
- **文件结构**：3 个生物重复 `.xls` 文件，每个都有 `Meta`、`Data Columns`、`Data` 三个 sheet，合并共 **`36,912` 条原始孔记录**（`12,420 + 12,072 + 12,420`）。
- **每行是什么**：一个 384 孔板上的单个孔，记录了 `CellLine`、`BRDid`（Broad 化合物 ID）、`LJP Library Plate`（`LJP5` 或 `LJP6`）、`Small Mol Concentration`（6 档浓度：`0.04, 0.12, 0.37, 1.11, 3.33, 10.0 uM`）、`Timepoint`（全部 `72 h`）以及成像核计数：
  - `Total Cell Count Before Treatment`（$N_0$，中位数 `1235.4`，范围 `513–2543`）
  - `Total Cell Count After Treatment`（$N_t$，中位数 `2692.0`，范围 `1–12883`）
  - `Total Control Cell Count`（$N_c$，中位数 `4156.2`，范围 `3125–11156`）
  - `Relative Cell Count`（$N_t / N_c$，中位数 `0.6162`，范围 `0.0002–1.3087`）
  - `Relative Growth Rate After Treatment`（$(N_t - N_0)/(N_c - N_0)$，中位数 `0.4733`，范围 `-0.7592–1.5096`；其中 **`7,668` 孔 / `20.77%` 为负值**，代表细胞数比给药前更少，即**净杀伤**）。
- **覆盖规模**：`6` 株乳腺细胞（`BT20, HS578T, MCF10A, MCF7, MDAMB231, SKBR3`）× `107` 个药物。
- **看一条真实原孔样本（`Screen20237`, `LJP6`, `BT20 + afatinib (BRD-A58767537) + 0.04 uM`）**：
  - 3 个技术重复孔：$N_0 = 1657.74$、$N_c = 5505.67$、处理后 $N_t = 3973 / 4132 / 3836$ → 相对细胞数分别为 `0.72162 / 0.75050 / 0.69674`（均值 `0.72295`）。

### 3.2 `LINCS_GSE70138`（Phase II L1000 表达数据）
- **基因表 `gene_info`**：共 `12,328` 个基因，其中仅 **`978` 个是实际探针测量的 landmark 基因（`pr_is_lm == 1`）**，其余 `11,350` 个是由线性模型推断出来的。因此我们始终只用这 **978 个实测基因**，避免把推断噪声当真值。
- **签名表 `sig_info`**：共 `118,050` 条汇总签名，其中小分子处理 `trt_cp` 为 `107,404` 条（`24 h`: `101,360`、`3 h`: `4,950`、`6 h`: `1,094`）；与 HMS 对应的 `LJP005/LJP006` 库板共 `14,832` 条（`24 h`: `9,882`、`3 h`: `4,950`）。
- **实例表 `inst_info`**：共 `345,976` 个单孔实例（`314,422` 个 `trt_cp`、`19,258` 个 `ctl_vehicle` DMSO 溶剂对照），记录了每条 Level 5 签名由哪些具体实验孔（`distil_id`）和哪块检测板（`det_plate`）汇总而来。
- **药物表 `pert_info`（发现陷阱 #1）**：`1,796` 个 `trt_cp` 化合物中有 `2` 个的 `canonical_smiles` 是 `restricted` 或缺失（进入 HMS 匹配时会有 1 个受限结构药无法生成分子指纹，必须显式剔除）。
- **两大核心 GCTX 矩阵（HDF5 格式）**：
  - `Level 3`（`15.90 GB`，shape `345,976 × 12,328`, `float32`）：单孔归一化表达值（抽样范围 `0.0–15.0`，中位数 `9.63`），用于逐实例核查有效性与同板 DMSO 对照数；
  - `Level 5`（`5.42 GB`，shape `118,050 × 12,328`, `float32`）：相对参考的标准化响应 z-score（抽样范围 `-10.0 ~ +10.0`，中位数 `~0.0`），正负表示上调/下调，这就是 `v0` 的 978 维表达监督标签。

### 3.3 `CCLE_DepMap_24Q4`（细胞基线 RNA-seq）
- **身份表 `Model.csv`**：`2,105` 条细胞模型记录，提供稳定的跨库身份键 `ModelID`（如 `ACH-000536`）、`RRID`、`OncotreeLineage`、`SangerModelID`。
- **基线表达矩阵 `OmicsExpressionProteinCodingGenesTPMLogp1.csv`**（`483.16 MB`）：**`1,673` 个有 RNA 的模型 × `19,193` 个蛋白编码基因**，列名格式为 `基因符号 (EntrezID)`（如 `TSPAN6 (7105)`），单位为 `log2(TPM + 1)`（抽样范围 `0.0–9.49`，中位数 `3.14`）。
- **发现陷阱 #2（为什么 `v0` 只有 5 株细胞而不是 6 株？）**：脚本直接逐一核对了 HMS 的 6 株细胞：

| HMS 细胞名 | DepMap `ModelID` | 谱系 / 疾病类型 | CCLE 24Q4 是否有 RNA 行 |
| :--- | :--- | :--- | :---: |
| `BT20` | `ACH-000536` | Breast / Invasive Breast Carcinoma | **是** |
| `HS578T` | `ACH-000148` | Breast / Invasive Breast Carcinoma | **是** |
| `MCF7` | `ACH-000019` | Breast / Invasive Breast Carcinoma | **是** |
| `MDAMB231` | `ACH-000768` | Breast / Invasive Breast Carcinoma | **是** |
| `SKBR3` | `ACH-000017` | Breast / Invasive Breast Carcinoma | **是** |
| **`MCF10A`** | `ACH-001357` | Breast / **Non-Cancerous**（非癌上皮） | **否（无 24Q4 RNA 行）** |

---

## 4. 代码、改动与如何手动运行

- **新增代码**：[inspect_raw_data.py](inspect_raw_data.py)（零外部状态依赖，只读扫描原始文件）
- **实测产出**：[raw_data_inventory.csv](raw_data_inventory.csv)、[raw_data_profile.json](raw_data_profile.json)
- **相对上一步改动**：本步为新分支 `step-by-step-repro` 的第 0 步起点，建立了迭代路线图 [iterations/README.md](../README.md) 与原始数据体检基准。
- **手动运行命令**：
  ```bash
  python iterations/00_raw_data_audit/inspect_raw_data.py
  # 如原始数据放在其他路径，可通过参数指定：
  python iterations/00_raw_data_audit/inspect_raw_data.py --raw-root /path/to/raw --unpacked-root /path/to/unpacked
  ```

---

## 5. 引出下一步（`iter-01`）

现在原始数据的字段、单位、两个核心陷阱（`MCF10A` 缺 CCLE RNA、受限 SMILES）以及库板/重复结构已经全部看清。下一步 [`iter-01`](../01_v0_preprocessing/README.md) 将动手写**最初版本（`v0`）的数据预处理代码**，把这 3 个原始数据源一步步洗成 `5 细胞 × 106 药 × 3,364 条观测` 的干净训练数组与严格防泄漏划分。
