# 数据说明

## 已入库（`data/processed/`）

| 文件 | 内容 | 行数 |
|---|---|---|
| `ghsr1a_binding_clean.csv` | 清洗后的结合数据（IC50 + Ki），建模主数据 | 1355 |
| `ghsr1a_functional_clean.csv` | 清洗后的功能数据（EC50），单独建模用 | 832 |

`*_clean.csv` 的主要列：`canonical_smiles`、`pAffinity` / `pEC50`、`MW`、`logP`、
`HBD`、`HBA`、`TPSA`、`QED`、`scaffold`、`peptide_like`、`rotatable_bonds`。

**结合数据与功能数据不合并建模** —— GHSR1a 组成性活性约 50%（GPCR 中最高之一），
功能测定的 IC50 表示"抑制基础活性"，与 EC50 的"激活"语义相反，混入同一模型会造成冲突。

## 未入库（`data/raw/`，可重新生成）

原始数据由 `notebooks/01_数据集构建.ipynb` 从 ChEMBL 下载，耗时约 3–8 分钟。
首次运行会自动缓存到 `data/raw/`，之后重跑直接读取。

```bash
cd notebooks && jupyter lab
# 运行 01_数据集构建.ipynb
```

## 数据来源

| 来源 | 用途 | 标识 |
|---|---|---|
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | 化合物活性数据 | 靶点 `CHEMBL4616`（人源 GHSR1a，UniProt Q92847） |
| [RCSB PDB](https://www.rcsb.org/) | 受体结构 | `9UY3`（2.52 Å，anamorelin 共晶） |
| [UniProt](https://www.uniprot.org/) | 序列与跨膜区标注 | `Q92847` |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | 配体结构校验 | anamorelin，C31H42N6O3 |
| [GEO](https://www.ncbi.nlm.nih.gov/geo/) | 转录组（探索性） | `GSE307669`、`GSE308625` |

## 重新获取受体结构

```bash
# 9UY3 只有 mmCIF 格式（结构过大，RCSB 不提供传统 PDB）
curl -O https://files.rcsb.org/download/9UY3.cif

# 然后提取受体链与共晶配体
python scripts/structure/prepare_9UY3.py
```

## 关于物种

活性数据已限定 `target_organism = Homo sapiens`。
不要混入 `CHEMBL3278`（大鼠）或 `CHEMBL3428`（小鼠）——跨物种活性不可合并。
