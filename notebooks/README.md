# Notebooks

按编号顺序运行。每个 notebook 的第一个代码单元会自动定位项目根目录，
因此在仓库内任何位置打开都能正确读写——文件会写入 `data/` 与 `results/` 对应目录。

| Notebook | 内容 | 耗时 |
|---|---|---|
| `01_数据集构建.ipynb` | ChEMBL 取数 → 清洗 → 探索 → 保存 | 首次 3–8 分钟（之后读缓存） |
| `02_数据切分对比实验.ipynb` | 六种切分方式对比 + 新颖度分层诊断 | 约 1 分钟 |
| `03_QSAR建模.ipynb` | 特征工程 → 模型网格 → 适用域判断 → 批量预测 | 约 3 分钟 |

## 前置条件

```bash
conda activate aidd
cd notebooks
jupyter lab
```

需要 `data/processed/ghsr1a_binding_clean.csv` 存在（已入库）。
若想从零重建，先跑 01。

## 运行提示

- **内核忙碌时右上角是实心圆** —— 在正常工作，不是在死机
- **中断** 点工具栏黑色方块 `■`，不要直接关窗口
- **关闭** 在终端按 `Control + C`，输入 `y` 回车

## 关于路径

第一个代码单元打印类似下面的内容，说明路径已正确初始化：

```
项目根目录: /path/to/GHSR1a-AIDD
  数据 → data/processed
  图表 → results/figures
```

若显示"未定位到项目根目录，使用当前工作目录（扁平模式）"，说明 notebook 被复制到了
仓库之外。此时它会在自身所在目录读写——把 `data/processed/` 下的 CSV 复制过去即可。
