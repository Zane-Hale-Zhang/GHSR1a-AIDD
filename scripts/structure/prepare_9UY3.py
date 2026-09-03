#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GHSR1a 结构准备 · 从 mmCIF 提取受体与共晶配体
================================================
不依赖任何第三方包，纯标准库实现（因为 9UY3 只有 mmCIF 格式，没有传统 PDB）。

功能：
  1. 下载 9UY3.cif（如已存在则跳过）
  2. 列出结构里所有链与配体，看清这个复合物到底由什么组成
  3. 提取受体链 R（去掉 G 蛋白、纳米抗体、胆固醇）
  4. 提取共晶配体 A1ESD（anamorelin）单独存盘
  5. 计算配体质心 —— 这就是对接盒子的中心
  6. 建议对接盒子尺寸

用法：
    python prepare_9UY3.py
"""

import os, sys, urllib.request
from collections import OrderedDict, defaultdict

OUTDIR = os.path.dirname(os.path.abspath(__file__))
CIF = os.path.join(OUTDIR, '9UY3.cif')
URL = 'https://files.rcsb.org/download/9UY3.cif'

RECEPTOR_CHAIN = 'R'      # GHSR1a 本体
LIGAND_COMP = 'A1ESD'     # anamorelin
STRIP_COMPS = {'CLR', 'HOH', 'OLC', 'SO4', 'PO4', 'GOL', 'NAG', 'PLM'}


def sep(t='-', n=74):
    print(t * n)


def download():
    if os.path.exists(CIF) and os.path.getsize(CIF) > 100000:
        print(f'  已存在 {os.path.basename(CIF)} '
              f'({os.path.getsize(CIF)/1024/1024:.1f} MB)，跳过下载')
        return
    print(f'  正在下载 {URL} ...')
    urllib.request.urlretrieve(URL, CIF)
    print(f'  完成 ({os.path.getsize(CIF)/1024/1024:.1f} MB)')


def parse_atom_site(path):
    """解析 mmCIF 的 _atom_site 循环，返回 dict 列表"""
    with open(path, encoding='utf-8') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == 'loop_':
            j = i + 1
            cols = []
            while j < len(lines) and lines[j].strip().startswith('_atom_site.'):
                cols.append(lines[j].strip().split('.', 1)[1])
                j += 1
            if 'group_PDB' in cols:
                data, k = [], j
                while k < len(lines):
                    s = lines[k].strip()
                    if s in ('#', 'loop_', '') or s.startswith('_'):
                        break
                    p = s.split()
                    if len(p) == len(cols):
                        data.append(dict(zip(cols, p)))
                    k += 1
                return data
            i = j
        else:
            i += 1
    raise RuntimeError('未找到 _atom_site 循环')


def fmt_pdb(a, serial):
    """按标准 PDB 固定列格式输出一行"""
    name = a.get('label_atom_id', 'C')
    elem = a.get('type_symbol', 'C')
    if len(name) < 4 and len(elem) == 1:
        name = ' ' + name
    name = name.ljust(4)[:4]
    alt = a.get('label_alt_id', '.')
    alt = ' ' if alt in ('.', '?') else alt
    resn = a.get('label_comp_id', 'UNK')[-3:].rjust(3)
    ch = a.get('auth_asym_id', 'A')
    try:
        resi = int(a.get('auth_seq_id', 1))
    except ValueError:
        resi = 1
    x, y, z = float(a['Cartn_x']), float(a['Cartn_y']), float(a['Cartn_z'])
    occ = float(a.get('occupancy', 1.0) or 1.0)
    bfac = float(a.get('B_iso_or_equiv', 0.0) or 0.0)
    return (f'ATOM  {serial:5d} {name}{alt:1s}{resn:>3s} {ch:1s}{resi:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}          {elem:>2s}')


def main():
    sep('=')
    print('  GHSR1a 结构准备 · 9UY3 (anamorelin 复合物, 2.52 Å)')
    sep('=')

    print('\n[1/6] 获取结构文件')
    download()

    print('\n[2/6] 解析 mmCIF')
    atoms = parse_atom_site(CIF)
    print(f'  共 {len(atoms)} 个原子记录')

    # ---------- 结构组成概览 ----------
    print('\n[3/6] 这个复合物由什么组成')
    chains = defaultdict(lambda: {'res': OrderedDict(), 'het': OrderedDict()})
    for a in atoms:
        ch = a.get('auth_asym_id', '?')
        if a.get('group_PDB') == 'HETATM':
            comp = a.get('label_comp_id', '?')
            chains[ch]['het'].setdefault(comp, 0)
            chains[ch]['het'][comp] += 1
        else:
            try:
                chains[ch]['res'][int(a.get('auth_seq_id', 0))] = \
                    a.get('label_comp_id', '?')
            except ValueError:
                pass

    print(f"\n  {'链':<4}{'残基数':>7}{'范围':>14}   说明")
    print('  ' + '-' * 60)
    NOTES = {'R': '★ GHSR1a 受体本体（对接目标）',
             'A': 'Gαq 亚基（信号复合物，移除）',
             'B': 'Gβ 亚基（移除）',
             'G': 'Gγ 亚基（移除）',
             'N': 'Nb35 纳米抗体（稳定用，移除）'}
    for ch in sorted(chains):
        info = chains[ch]
        if not info['res']:
            continue
        rs = sorted(info['res'])
        note = NOTES.get(ch, '')
        print(f"  {ch:<4}{len(rs):>7}{f'{rs[0]}–{rs[-1]}':>14}   {note}")
        for comp, cnt in info['het'].items():
            tag = '（胆固醇，移除）' if comp == 'CLR' else '（共晶配体 ★保留）'
            print(f"         └ HETATM {comp:<6} {cnt:3d} 原子 {tag}")

    # ---------- 提取受体 ----------
    print(f'\n[4/6] 提取受体链 {RECEPTOR_CHAIN}')
    receptor = [a for a in atoms
                if a.get('auth_asym_id') == RECEPTOR_CHAIN
                and a.get('group_PDB') == 'ATOM']
    out_r = os.path.join(OUTDIR, '9UY3_receptor.pdb')
    with open(out_r, 'w') as f:
        f.write('REMARK   GHSR1a receptor chain R extracted from 9UY3\n')
        f.write('REMARK   2.52 A cryo-EM, anamorelin complex, G protein removed\n')
        for i, a in enumerate(receptor, 1):
            f.write(fmt_pdb(a, i) + '\n')
        f.write('END\n')
    nres = len({int(a['auth_seq_id']) for a in receptor})
    print(f'  → {out_r}')
    print(f'    {len(receptor)} 原子 / {nres} 残基')

    # ---------- 提取配体 ----------
    print(f'\n[5/6] 提取共晶配体 {LIGAND_COMP}')
    lig = [a for a in atoms
           if a.get('label_comp_id') == LIGAND_COMP]
    if not lig:
        print(f'  ✗ 未找到 {LIGAND_COMP}，请检查化学组分 ID')
        return
    out_l = os.path.join(OUTDIR, '9UY3_ligand_A1ESD.pdb')
    with open(out_l, 'w') as f:
        f.write(f'REMARK   anamorelin ({LIGAND_COMP}) from 9UY3\n')
        f.write('REMARK   用途: 定义口袋中心 + redocking 验证\n')
        for i, a in enumerate(lig, 1):
            f.write(fmt_pdb(a, i) + '\n')
        f.write('END\n')
    print(f'  → {out_l}')
    print(f'    {len(lig)} 原子')

    # ---------- 口袋与盒子 ----------
    print('\n[6/6] 对接盒子参数')
    xs = [float(a['Cartn_x']) for a in lig]
    ys = [float(a['Cartn_y']) for a in lig]
    zs = [float(a['Cartn_z']) for a in lig]
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    ex, ey, ez = max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)
    suggest = max(18.0, round(max(ex, ey, ez) + 10))

    print(f'\n  共晶配体 anamorelin 的尺寸: {ex:.1f} × {ey:.1f} × {ez:.1f} Å')
    print(f'\n  ┌─── 对接盒子（可直接用于 AutoDock Vina）─────────────')
    print(f'  │ center_x = {cx:.3f}')
    print(f'  │ center_y = {cy:.3f}')
    print(f'  │ center_z = {cz:.3f}')
    print(f'  │ size_x   = {suggest}')
    print(f'  │ size_y   = {suggest}')
    print(f'  │ size_z   = {suggest}')
    print(f'  └───────────────────────────────────────────────')
    print(f'\n  说明: 盒子中心取共晶配体的质心，这是最可靠的口袋定义方式。')
    print(f'        尺寸 = 配体最长边 ({max(ex,ey,ez):.1f} Å) + 10 Å 余量，')
    print(f'        保证配体在盒内有足够旋转空间。')

    sep('=')
    print('  完成。下一步:')
    print('    1. 用 PyMOL 打开 9UY3.cif，核对上面的判断')
    print('    2. 对受体加氢、确定质子化状态（pdbfixer 或 H++）')
    print('    3. 把 A1ESD 重新对接回去（redocking），验证流程正确')
    print('       —— RMSD < 2 Å 才算合格，这是整个对接工作流的校准步骤')
    sep('=')


if __name__ == '__main__':
    main()
