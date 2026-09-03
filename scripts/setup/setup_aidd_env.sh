#!/bin/bash
# ==========================================================
# AIDD 起步环境安装脚本 (macOS Apple Silicon 专用)
# 适用: Apple M 系列芯片 + macOS
# 作用: 安装 Miniconda, 创建 aidd 环境, 装好第一阶段所需全部库
# 可重复运行, 已存在的内容会自动跳过
# ==========================================================

set -e

CONDA_DIR="$HOME/miniconda3"
ENV_NAME="aidd"
LOG="$HOME/aidd_setup.log"

say() { echo ""; echo "──────────────────────────────────────"; echo "  $1"; echo "──────────────────────────────────────"; }

exec > >(tee -a "$LOG") 2>&1

say "步骤 1/6 · 下载 Miniconda (Apple Silicon 原生版)"
if [ -f "$CONDA_DIR/bin/conda" ]; then
    echo "Miniconda 已存在, 跳过下载"
else
    echo "正在下载, 约 100MB, 请稍候..."
    curl -L --max-time 900 -o /tmp/miniconda.sh \
        https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
    echo "下载完成, 开始安装到 $CONDA_DIR"
    bash /tmp/miniconda.sh -b -p "$CONDA_DIR"
    rm -f /tmp/miniconda.sh
    echo "Miniconda 安装完成"
fi

say "步骤 2/6 · 初始化 shell"
"$CONDA_DIR/bin/conda" init zsh >/dev/null 2>&1
echo "已写入 zsh 配置"

say "步骤 3/6 · 配置软件源 (纯 conda-forge, 完全绕开 Anaconda 商业频道)"
"$CONDA_DIR/bin/conda" config --remove channels defaults 2>/dev/null || true
"$CONDA_DIR/bin/conda" config --remove-key default_channels 2>/dev/null || true
"$CONDA_DIR/bin/conda" config --add channels conda-forge
"$CONDA_DIR/bin/conda" config --add default_channels https://conda.anaconda.org/conda-forge
"$CONDA_DIR/bin/conda" config --set channel_priority strict
"$CONDA_DIR/bin/conda" config --set auto_activate_base false
echo "软件源配置完成: conda-forge, strict 优先级"
echo ""
echo "  【重要说明】新版 conda 使用 Anaconda 商业频道前必须先接受其服务条款,"
echo "  否则建环境时会报错: CondaToSNonInteractiveError。"
echo "  本脚本把 default_channels 也一并改到 conda-forge, 彻底绕开该问题 ——"
echo "  既免去了条款授权的法律灰度, 也避免了 defaults 与 conda-forge 混用"
echo "  导致的底层库冲突(这是 RDKit 崩溃的头号原因)。"

say "步骤 4/6 · 创建 $ENV_NAME 环境 (Python 3.10)"
if "$CONDA_DIR/bin/conda" env list | grep -qE "^${ENV_NAME}\s"; then
    echo "环境 $ENV_NAME 已存在, 跳过创建"
else
    "$CONDA_DIR/bin/conda" create -n "$ENV_NAME" python=3.10 -y
    echo "环境创建完成"
fi

say "步骤 5/6 · 安装核心库"
rm -f "$CONDA_DIR/envs/$ENV_NAME"/.tmp.* 2>/dev/null || true
echo "包含: RDKit / pandas / numpy / scikit-learn / matplotlib / seaborn / jupyterlab"
echo "首次安装需要下载较多文件, 约 3-10 分钟, 请耐心等待..."
"$CONDA_DIR/bin/conda" install -n "$ENV_NAME" -y \
    rdkit \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    jupyterlab \
    tqdm \
    requests

say "步骤 6/6 · 验证安装"
"$CONDA_DIR/envs/$ENV_NAME/bin/python" - <<'PYEOF'
print("")
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED
import pandas, sklearn, numpy

smiles = 'CC(=O)Oc1ccccc1C(=O)O'   # 阿司匹林
mol = Chem.MolFromSmiles(smiles)
print("  测试分子 (阿司匹林):", smiles)
print("  原子数            :", mol.GetNumAtoms())
print("  分子量            :", round(Descriptors.MolWt(mol), 2))
print("  logP              :", round(Crippen.MolLogP(mol), 2))
print("  QED 类药性评分    :", round(QED.qed(mol), 3))
print("")
import rdkit, sys
print("  Python 版本 :", sys.version.split()[0])
print("  RDKit 版本  :", rdkit.__version__)
print("  pandas 版本 :", pandas.__version__)
print("  sklearn 版本:", sklearn.__version__)
print("")
print("  全部通过, 环境可用")
PYEOF

say "设置终端默认进入 $ENV_NAME 环境"
if ! grep -q "conda activate $ENV_NAME" "$HOME/.zshrc" 2>/dev/null; then
    echo "" >> "$HOME/.zshrc"
    echo "# AIDD 学习环境: 打开终端自动进入" >> "$HOME/.zshrc"
    echo "conda activate $ENV_NAME" >> "$HOME/.zshrc"
    echo "已设置"
else
    echo "已设置过, 跳过"
fi

say "安装全部完成"
echo ""
echo "  ⚠️  必须做: 完全退出终端再重新打开 (Command + Q), 或执行:"
echo "        source ~/.zshrc"
echo ""
echo "  原因: 本脚本刚刚才把 conda 的初始化代码写进 ~/.zshrc,"
echo "        而你当前这个终端窗口是在写入之前打开的, 它对此一无所知。"
echo "        不重开的话, 敲任何 conda / jupyter 命令都会报 command not found。"
echo ""
echo "  成功标志: 命令行最左边出现 (aidd) 字样"
echo "  安装日志: $LOG"
echo ""
