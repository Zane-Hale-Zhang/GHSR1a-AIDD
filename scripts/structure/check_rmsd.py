#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD 检查 · 对接结果 vs 晶体构象（支持多构象批量比较）
================================================================
用法:
    python check_rmsd.py 参考.pdb 结果1.pdb [结果2.pdb ...]
    python check_rmsd.py 参考.pdb pose*.pdb

为什么支持多个:
    对接软件输出的构象按「打分」排序，但打分最好的不一定是位置最准的。
    标准做法报告两个指标:
      · Top-1  RMSD  —— 打分第一的构象准不准（衡量打分函数）
      · Best-of-N    —— N 个里最接近晶体的（衡量采样能力）
    若 Best-of-N 合格而 Top-1 不合格，说明采样够、打分排序不准，
    这种对接仍可用于筛选，只是不能盲信排名。

判定标准:
    RMSD < 2.0 Å   合格
    RMSD 2.0-3.0 Å 勉强
    RMSD > 3.0 Å   失败
"""

import sys, os, glob
import numpy as np

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


def sep(t='-', n=66):
    print(t * n)


def load_heavy_coords(path):
    """读重原子坐标（原子序数 > 1）。
    注意: 对接软件输出的 PDB 常把键序全标成单键（CONECT 不含键级），
    此时 RDKit 的 removeHs 会失效，故这里按原子序数手动过滤氢。
    """
    if HAS_RDKIT:
        ext = os.path.splitext(path)[1].lower()
        try:
            if ext in ('.sdf', '.sd'):
                m = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
            elif ext == '.mol2':
                m = Chem.MolFromMol2File(path, removeHs=False, sanitize=False)
            else:
                m = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
            if m is not None and m.GetNumAtoms() > 0 and m.GetNumConformers() > 0:
                conf = m.GetConformer()
                idx = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
                xyz = np.array([list(conf.GetAtomPosition(i)) for i in idx])
                return xyz, m
        except Exception:
            pass

    # 回退：PDB 固定列位，按元素列过滤氢
    xyz = []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            elem = line[76:78].strip().upper()
            if elem == 'H':
                continue
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz), None


def distmat(X):
    return np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)


def order_consistency(ref_xyz, dock_xyz):
    """用内部距离矩阵相关性检查原子顺序是否一致"""
    if len(ref_xyz) != len(dock_xyz):
        return None
    Dr, Dd = distmat(ref_xyz), distmat(dock_xyz)
    iu = np.triu_indices(len(ref_xyz), k=1)
    return float(np.corrcoef(Dr[iu], Dd[iu])[0, 1])


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    ref_path = sys.argv[1]
    dock_paths = []
    for arg in sys.argv[2:]:
        if '*' in arg or '?' in arg:
            dock_paths.extend(sorted(glob.glob(arg)))
        else:
            dock_paths.append(arg)
    dock_paths = [p for p in dock_paths if os.path.abspath(p) != os.path.abspath(ref_path)]

    if not dock_paths:
        print('✗ 没有可比较的结果文件')
        sys.exit(1)

    if not os.path.exists(ref_path):
        print(f'✗ 参考文件不存在: {ref_path}')
        sys.exit(1)

    sep('=')
    print('  RMSD 检查 · 对接结果 vs 晶体构象')
    sep('=')
    print(f'\n  参考 (晶体): {os.path.basename(ref_path)}')
    print(f'  待比较     : {len(dock_paths)} 个构象')

    ref_xyz, _ = load_heavy_coords(ref_path)
    print(f'  参考重原子 : {len(ref_xyz)} 个')

    results = []
    for p in dock_paths:
        if not os.path.exists(p):
            print(f'  ⚠ 跳过（不存在）: {p}')
            continue
        d_xyz, _ = load_heavy_coords(p)
        if len(d_xyz) != len(ref_xyz):
            print(f'  ⚠ 跳过（原子数 {len(d_xyz)} ≠ {len(ref_xyz)}）: {os.path.basename(p)}')
            continue
        rms = float(np.sqrt(np.mean(np.sum((ref_xyz - d_xyz) ** 2, axis=1))))
        corr = order_consistency(ref_xyz, d_xyz)
        results.append((os.path.basename(p), rms, corr))

    if not results:
        print('\n✗ 没有可比较的构象')
        sys.exit(1)

    # 汇总表
    print('\n' + '-' * 66)
    print(f"  {'构象文件':<34}{'RMSD':>9}{'顺序检查':>11}")
    print('  ' + '-' * 62)
    for name, rms, corr in results:
        flag = '✓' if rms < 2.0 else ('~' if rms < 3.0 else '✗')
        cs = f'{corr:.3f}' if corr is not None else '  -  '
        warn = '' if (corr is None or corr > 0.85) else '  ⚠顺序存疑'
        print(f'  {name[:33]:<34}{rms:>7.2f} Å{cs:>10}{warn}')
    print('  ' + '-' * 62)

    best_name, best_rms, _ = min(results, key=lambda x: x[1])
    top1_rms = results[0][1]

    print()
    sep('=')
    print('  结果判定')
    sep('=')
    print(f'\n  Top-1（打分最好的构象）  : {top1_rms:.2f} Å  ← 衡量打分函数准不准')
    print(f'  Best-of-{len(results)}（最接近晶体的）  : {best_rms:.2f} Å  ← 衡量采样能力够不够')
    print(f'    （最佳构象: {best_name}）')

    print()
    if best_rms < 2.0 and top1_rms < 2.0:
        print('  ✓ 完全合格 —— 打分和采样都对，可放心用于虚拟筛选')
    elif best_rms < 2.0:
        print('  ⚠ 采样合格、打分排序不准')
        print('    含义: 对接能找到正确的结合模式，但排不进第一名。')
        print('    用法: 可用于筛选，但别只看排名第一的构象，')
        print('          建议取前几个构象人工检查或做二次打分（rescoring）。')
    elif best_rms < 3.0:
        print('  ~ 勉强 —— N 个构象里最好的仍有偏差')
        print('    建议: 增加构象数、提高 exhaustiveness，或检查口袋是否选对')
    else:
        print('  ✗ 失败 —— 采样和打分都有问题')
        print('    排查: 口袋是否选对 → 配体结构 → 受体质子化 → 分子柔性是否过大')

    avg_corr = np.mean([r[2] for r in results if r[2] is not None])
    if avg_corr < 0.85:
        print()
        print(f'  ⚠ 注意：原子顺序一致性偏低 (r={avg_corr:.2f})，RMSD 可能不可靠。')
        print('    常见原因：对接软件重排了原子顺序，或氢的处理方式不同。')
    sep('=')


if __name__ == '__main__':
    main()
