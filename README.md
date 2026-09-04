# GHSR1a 药物发现数据科学项目

> 从公共数据到可验证预测模型的完整计算药物发现流水线——
> 包含数据集构建、受体结构解析、对接流程验证、QSAR 建模与方法学严谨性检验。

**靶点**：GHSR1a（生长激素促分泌素受体 1a，ghrelin 受体），A 类 GPCR
**结构模板**：9UY3（2.52 Å，anamorelin 共晶）
**数据来源**：ChEMBL `CHEMBL4616` · RCSB PDB · UniProt `Q92847`

---

## 技术栈

`Python 3.10` · `RDKit` · `ChEMBL API` · `PyMOL` · `PDBFixer` · `OpenBabel` · `gemmi`
· `AutoDock Vina`（CB-Dock2）· `scikit-learn` · `XGBoost` · `LightGBM` · `pandas` · `matplotlib`

---

## 项目亮点

**完整链路** — 不是孤立的模型，而是从公共数据库取数到可部署预测函数的端到端流程。

**方法学严谨** — 报告 scaffold split 而非随机切分；多随机种子取平均；
逐 seed 对比以区分真实改善与噪声；装包前 dry-run 查依赖冲突。

**对自身局限有清晰界定** — 通过实验确定了模型的失效临界点
（测试分子与训练集相似度 < 0.6 时失效，域外 R² 为负），
并将其固化为预测函数的适用域判断，而非只报告一个漂亮的分数。

**包含诚实的负面结果** — Redocking 未达标（RMSD 3.30 Å）、
XGBoost/LightGBM 未能超越基线模型，均完整记录并给出诊断。

---

## 快速开始

```bash
# 1. 环境（macOS Apple Silicon；其他平台见 environment.yml 顶部注释）
conda env create -f environment.yml
conda activate aidd

# 若已装过 conda，只需补频道配置：
#   conda config --remove-key default_channels
#   conda config --add default_channels https://conda.anaconda.org/conda-forge
#   conda config --set channel_priority strict
```

```bash
# 2. 数据（清洗后数据已入库，可直接开始）
ls data/processed/
#   ghsr1a_binding_clean.csv      1355 个化合物
#   ghsr1a_functional_clean.csv    832 个化合物

# 想从零重建，运行：notebooks/01_数据集构建.ipynb（首次下载 3–8 分钟）
```

```bash
# 3. 结构准备
cd scripts/structure
python prepare_9UY3.py        # 纯标准库 mmCIF 解析，提取受体与共晶配体，算对接盒子
python analyze_pocket.py      # 口袋残基与跨膜区归属
python add_hydrogens.py       # 加氢 + 质子化状态检查 + 盐桥复测
```

```bash
# 4. 建模与评估
cd ../modeling
python modeling_grid.py       # 特征 × 模型网格对比
python boosting_comparison.py # XGBoost / LightGBM 对比

# 或运行 notebooks/03_QSAR建模.ipynb 查看完整分析与可视化
```

---

## 项目结构

```
GHSR1a-AIDD/
├── README.md                        本文件
├── LICENSE                          MIT
├── environment.yml                  可复现环境（含踩坑注释）
├── requirements.txt                 pip 备选
│
├── notebooks/                       分析主线，按编号运行
│   ├── 01_数据集构建.ipynb            ChEMBL 取数 → 清洗 → 探索 → 保存
│   ├── 02_数据切分对比实验.ipynb      六种切分方式 + 新颖度分层诊断
│   └── 03_QSAR建模.ipynb             特征工程 → 模型网格 → 适用域 → 批量预测
│
├── scripts/
│   ├── structure/                   结构生物学
│   │   ├── prepare_9UY3.py            mmCIF 解析器 + 受体/配体提取（无第三方依赖）
│   │   ├── analyze_pocket.py          口袋残基与 TM/ECL 归属
│   │   ├── add_hydrogens.py           加氢 + 质子化 + 盐桥复测
│   │   ├── check_rmsd.py              多构象 RMSD + 原子顺序自检
│   │   └── view_9UY3.pml              PyMOL 可视化
│   ├── modeling/
│   │   ├── modeling_grid.py           特征 × 模型网格
│   │   ├── boosting_comparison.py     XGBoost / LightGBM 对比
│   │   └── split_comparison_v3.py     切分严格度光谱
│   └── setup/                        环境安装脚本（幂等）
│
├── data/
│   ├── processed/                    清洗后数据（入库）
│   └── raw/                          原始数据（不入库，可重新下载）
│
├── results/
│   ├── figures/                      6 张分析图
│   ├── tables/                       网格结果 + 口袋残基清单
│   └── structures/                   加氢后受体、配体、redocking 结果
│
└── docs/
    ├── index.html                    项目作品集（Pages 入口，自包含可直接打开/打印）
    ├── AMPK机制最小实验方案.md        两周闭环的机制设计（含决策树与应急方案）
    ├── 四周执行表_AMPK机制.md         按天排的执行表，含 Go/No-Go 决策点
    ├── 试剂采购清单.md                按到位时间分组，含已核实货号
    ├── protocol_表面CD36检测.md       表面 CD36 定量完整 protocol
    ├── 方法学备忘录.md                 全部踩坑记录与解决方案
    └── 组学探索笔记.md                 转录组方向的探索与证伪
```

---

## 主要结果

### 数据集

| 项目 | 数值 |
|---|---|
| ChEMBL 靶点 | CHEMBL4616（人源 GHSR1a） |
| 原始活性记录 | 3455 条（IC50 1440 / Ki 755 / EC50 1259） |
| 清洗后结合数据 | 1355 个（小分子 889 / 肽类 466） |
| 清洗后功能数据 | 832 个 |
| 骨架多样性 | 595 种骨架，2.3 化合物/骨架 |
| SAR 强度（骨架内一致性） | 0.687（强） |
| 活性悬崖比例 | 8.8%（低，数据干净） |

### 结构解析（9UY3，2.52 Å）

| 项目 | 结果 |
|---|---|
| 受体链 R | 283 残基（37–334） |
| 口袋残基（配体 5 Å 内） | 27 个（跨膜 18 + 胞外环 8 + N 端 1） |
| 最近接触残基 | Glu124³·³³（2.76 Å） |
| E124–R283 盐桥 | 2.56 Å（加氢后复测不变） |
| 配体校验 | C31N6O3 = 40 重原子，与 PubChem anamorelin 吻合 |

**分叉口袋**：E124–R283 盐桥将口袋劈为 Cavity I（肽链结合）与
Cavity II（疏水，ghrelin 辛酰基伸入）——这是 GHSR1a 在已解析 GPCR 中独有的结构特征，
也从结构层面解释了"为什么 ghrelin 必须酰化才能激活受体"。

自写脚本算出的 6 个最近残基全部是文献中突变实验证实过的关键残基。

### QSAR 建模（scaffold split，5 seeds）

| 特征 | 模型 | R² | MAE | Spearman |
|---|---|---|---|---|
| Morgan + 描述符 | GradientBoosting | **0.453** | 0.599 | **0.665** |
| Morgan + 描述符 | XGBoost（调参） | 0.448 | **0.582** | 0.647 |
| Morgan | RandomForest | 0.401 | 0.614 | 0.616 |
| Morgan | LightGBM（默认） | 0.323 | 0.634 | 0.574 |
| Morgan | **ExtraTrees** | **0.061** | 0.773 | 0.492 |

基线（永远预测均值）MAE = 0.887 → 最佳模型 0.599，**误差下降 32%**。

ExtraTrees 崩盘的原因：它在每个分裂点随机选阈值，而 Morgan 指纹是 2048 维
稀疏二值向量，随机阈值在这种特征上几乎提不出信息。
**模型没有万能的，得看特征的脾气来选。**

### 适用域判断

| 组别 | 占比 | MAE | R² |
|---|---|---|---|
| 适用域内（相似度 ≥ 0.6） | 86% | 0.605 | 0.444 |
| 适用域外（相似度 < 0.6） | 14% | 0.968 | **−0.655** |

**域外 R² 为负：在那里使用模型不如直接猜测平均值。**

这个临界点换更强的模型也改不了——它是数据的性质，不是模型的问题。

### 特征重要性（排列重要性，描述符部分）

`BertzCT 0.053` > **`logP 0.045`** > `MW 0.023` > `FracCsp3 0.021` > 其余 ≈ 0

logP 排第二，而 Cavity II 恰好是疏水腔（ghrelin 辛酰基伸入处）。
模型从未被告知任何口袋结构信息，却从 889 个分子的活性数据中
"重新发现"了亲脂性的关键作用——**统计模型与结构生物学两条独立路径互证。**

---

## 关键设计决策

**结合数据与功能数据分离** — GHSR1a 组成性活性约 50%（GPCR 中最高之一），
功能测定中的 IC50 表示"抑制基础活性"，与 EC50 的"激活"语义相反，
混入同一模型会造成冲突。

**使用 scaffold split** — 随机切分会让同骨架分子同时出现在训练集与测试集，
模型只需识别骨架即可获得高分，性能虚高。
本项目实测随机切分 R² 虚高 0.071。

**先量化数据质量再建模** — SAR 强度 0.687 证明构效关系真实存在。
若该值接近 1.0，说明骨架不携带信息，任何模型都只是在拟合噪声。

**Redocking 用不叠合的 RMSD** — `GetBestRMS` 会先做最优叠合，
把整体平移消掉，得到的是"构象像不像"而非"有没有放回原位"。
实测整体平移 1 Å：逐点比较 1.73 Å，GetBestRMS 0.00 Å。

**装包前 dry-run 检查依赖** — 实际拦截过一次风险：`boost-cpp`/`swig`
会将 RDKit 从 2026.03.5 降级到 2024.03.6，`mdtraj` 会降级 numpy/pandas。
两者均已排除，编译型依赖一律放入独立环境。

**交叉验证关键结论** — 口袋残基用自写脚本算出 27 个，PyMOL 独立验证同为 27 个；
失效临界点 0.6 由 1780 个测试分子的分层诊断与最终模型的独立实测共同印证。

---

## 已知局限

- **样本量偏小**：889 个小分子，限制了模型复杂度与外推能力
- **数据噪声**：ChEMBL 为多实验室数据汇总，系统性误差无法通过算法消除
- **柔性分子对接**：anamorelin 约 10 个可旋转键，redocking RMSD 3.30 Å。
  口袋定位正确（质心差 0.78 Å），但柔性脂肪链朝向不准——
  打分可用于排序，具体结合模式不可全信
- **换模型已到头**：XGBoost / LightGBM 均未实质超越 GradientBoosting，
  进一步提升需依赖数据与特征，而非算法
- **GHSR 在心脏 bulk 转录组中几乎测不到**（FPKM 0.001，0/16 检出），
  需单细胞测序定位其细胞来源（详见 `docs/组学探索笔记.md`）

---

## 文档

| 文件 | 内容 |
|---|---|
| `docs/index.html` | 项目作品集，自包含无外部依赖，可直接打开或打印为 PDF |
| `docs/AMPK机制最小实验方案.md` | 两周闭环的 AMPK–CD36 机制设计（含决策树与应急方案） |
| `docs/四周执行表_AMPK机制.md` | 按天排的实验执行表，含 Go/No-Go 决策点 |
| `docs/试剂采购清单.md` | 按到位时间分组的试剂清单，含已核实货号 |
| `docs/protocol_表面CD36检测.md` | 表面 CD36 定量完整 protocol |
| `docs/方法学备忘录.md` | 全部踩坑记录：环境冲突、API 变化、RMSD 陷阱、模型选择 |
| `docs/组学探索笔记.md` | 转录组方向的调研、发现与证伪过程 |

---

## 数据来源

- [ChEMBL](https://www.ebi.ac.uk/chembl/) — 化合物活性数据
- [RCSB PDB](https://www.rcsb.org/) — 9UY3 等受体结构
- [UniProt](https://www.uniprot.org/) — Q92847 序列与跨膜区标注
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) — anamorelin 结构校验
- [GEO](https://www.ncbi.nlm.nih.gov/geo/) — HFpEF 心脏转录组（探索性分析）
- [CB-Dock2](https://cadd.labshare.cn/cb-dock2/php/index.php) — 在线分子对接

---

## 许可

MIT
