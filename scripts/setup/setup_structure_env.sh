#!/bin/bash
# ==========================================================
# 结构生物学工具安装脚本 (macOS Apple Silicon)
# 修订版 v2 —— 已排除会破坏主环境的包
#
# v1 的教训: 一次性装 8 个包导致依赖冲突, 且 libmamba solver
#            报错解释功能崩溃(bad_variant_access), 掩盖了真实原因。
#            逐个 dry-run 排查后发现两个"危险品":
#              · boost-cpp + swig → 把 RDKit 从 2026.03.5 降级到 2024.03.6
#              · mdtraj           → 把 numpy 2.2.6 降到 1.26.4、pandas 降到 2.2.2
#            这两个都已从主环境移除。
# ==========================================================

set -e

CONDA_DIR="$HOME/miniconda3"
ENV_NAME="aidd"
LOG="$HOME/aidd_structure_setup.log"

say() { echo ""; echo "──────────────────────────────────────"; echo "  $1"; echo "──────────────────────────────────────"; }

exec > >(tee -a "$LOG") 2>&1

say "安装结构生物学工具（已验证的安全组合）"
echo "  PyMOL        看结构、出图              ✓ 无降级"
echo "  PDBFixer     补缺失残基、加氢          ✓ 无降级"
echo "  OpenBabel    分子格式转换              ✓ 无降级"
echo "  gemmi        mmCIF ↔ PDB 转换          ✓ 无降级"
echo ""
echo "  刻意不装（会破坏主环境）:"
echo "  ✗ mdtraj       → 会把 numpy 2.2.6 降到 1.26.4、pandas 降到 2.2.2"
echo "  ✗ boost-cpp    → 会把 RDKit 2026.03.5 降到 2024.03.6"
echo "  ✗ swig         → 同上"
echo "                   这两个是给 Vina 编译用的，得不偿失。"
echo "                   需要 MDTraj 时单独建环境，不要污染 aidd。"
echo ""
echo "  首次安装约 3-8 分钟（PDBFixer 会带一个较大的 OpenMM）"

"$CONDA_DIR/bin/conda" install -n "$ENV_NAME" -y -c conda-forge \
    pymol-open-source \
    pdbfixer \
    openbabel \
    gemmi

say "验证安装"
"$CONDA_DIR/envs/$ENV_NAME/bin/python" - <<'PYEOF'
import importlib
for m, use in [('pymol', '结构可视化'),
               ('pdbfixer', '结构修补与加氢'),
               ('openbabel', '格式转换'),
               ('gemmi', 'mmCIF 处理'),
               ('rdkit', '化学信息学（应仍为 2026.03.5）')]:
    try:
        mod = importlib.import_module(m)
        ver = getattr(mod, '__version__', '')
        print(f"  ✓ {m:12s} {ver:14s} {use}")
    except Exception:
        print(f"  ✗ {m:12s} 未安装")
print()
import numpy, pandas
print(f"  numpy  {numpy.__version__}   (应保持 2.2.6)")
print(f"  pandas {pandas.__version__}  (应保持 2.3.3)")
PYEOF

say "命令行工具检查"
echo -n "  obabel : "; "$CONDA_DIR/envs/$ENV_NAME/bin/obabel" -V 2>&1 | head -1 || echo "不可用"
echo -n "  gemmi  : "; "$CONDA_DIR/envs/$ENV_NAME/bin/gemmi" --version 2>&1 | head -1 || echo "不可用"

say "可选 · 独立环境试装 AutoDock Vina"
echo "  说明: Vina 需要 Boost + SWIG 编译，这两个包会污染 aidd 环境。"
echo "        正确做法是建一个独立环境，失败也无所谓。"
echo ""
read -t 20 -p "  是否尝试？(20 秒后自动跳过) [y/N] " ans || ans="n"
if [[ "$ans" =~ ^[Yy]$ ]]; then
    "$CONDA_DIR/bin/conda" create -n docking python=3.10 -y 2>&1 | tail -3
    if "$CONDA_DIR/bin/conda" install -n docking -y -c conda-forge boost-cpp swig 2>&1 | tail -3; then
        if "$CONDA_DIR/envs/docking/bin/pip" install --no-cache-dir vina 2>&1 | tail -3; then
            echo "  ✓ Vina 装在 docking 环境，用法: conda activate docking"
        else
            echo "  ✗ Vina 编译失败 —— 用 Google Colab 跑对接即可"
        fi
    else
        echo "  ✗ 依赖装不上 —— 用 Google Colab 跑对接即可"
    fi
else
    echo "  已跳过。对接用 Google Colab 即可，本地只做结构准备。"
fi

say "完成"
echo ""
echo "  重开终端（或 source ~/.zshrc）后验证:"
echo "      pymol ~/aidd/structures/9UY3.cif"
echo ""
echo "  日志: $LOG"
echo ""
