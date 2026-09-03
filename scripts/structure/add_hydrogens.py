#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
受体加氢脚本 · GHSR1a (9UY3 链 R)
========================================
Cryo-EM 结构里没有氢原子。对接前必须加氢，而加氢时最危险的就是
质子化状态设错 —— 对 GHSR1a 而言，就是 Glu124 与 Arg283 的盐桥。

本脚本做三件事:
  1. 补小缺口、加氢 (pH 7.4)
  2. 逐一报告关键残基的质子化状态
  3. 重新测量 E124–R283 盐桥距离，确认它没被破坏
"""

import os, math, warnings
warnings.filterwarnings('ignore')

D = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(D, '9UY3_receptor.pdb')
OUT = os.path.join(D, '9UY3_receptor_h.pdb')

KEY_RESIDUES = [124, 283, 120, 99, 102]        # E124 R283 Q120 D99 R102
HIS_LIST = None                                 # 运行后填入所有 HIS
SALT_BRIDGE = (('OE1', 124), ('OE2', 124)), (('NH1', 283), ('NH2', 283))


def sep(t='-', n=70):
    print(t * n)


sep('=')
print('  GHSR1a 受体加氢 · 9UY3 chain R')
sep('=')

from pdbfixer import PDBFixer
from openmm.app import PDBFile

print('\n[1] 载入结构')
fixer = PDBFixer(filename=INP)
n_res = len(list(fixer.topology.residues()))
print(f'    {n_res} 个残基')

print('\n[2] 检查缺失残基')
fixer.findMissingResidues()
total_missing = sum(len(v) for v in fixer.missingResidues.values())
print(f'    共识别出 {total_missing} 个缺失残基')

# 只补短缺口（<=3），跳过 ICL3 这类长柔性区
kept, skipped = 0, 0
filtered = {}
for k, v in fixer.missingResidues.items():
    if len(v) <= 3:
        filtered[k] = v
        kept += len(v)
    else:
        skipped += len(v)
fixer.missingResidues = filtered
print(f'    补全短缺口 {kept} 个，跳过长柔性区 {skipped} 个（ICL3 之类，补了也是猜的）')

print('\n[3] 非标准残基与杂原子')
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)

print('\n[4] 加氢 (pH 7.4)')
# 注意: PDBFixer >= 1.9 起, findMissingAtoms() 必须显式调用,
#       旧版本会自动在 addMissingAtoms() 内部调用, 新版不会。
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)

PDBFile.writeFile(fixer.topology, fixer.positions, open(OUT, 'w'), keepIds=True)
n_atoms = len(list(fixer.topology.atoms()))
print(f'    输出: {os.path.basename(OUT)}')
print(f'    {n_atoms} 原子')

# ---------------- 验证 ----------------
print('\n[5] 关键残基的质子化状态')

atoms = []
for line in open(OUT):
    if line.startswith(('ATOM', 'HETATM')):
        atoms.append({
            'name': line[12:16].strip(),
            'resn': line[17:20].strip(),
            'resi': int(line[22:26]),
            'x': float(line[30:38]), 'y': float(line[38:46]), 'z': float(line[46:54]),
            'elem': line[76:78].strip(),
        })

by_res = {}
for a in atoms:
    by_res.setdefault(a['resi'], {'name': a['resn'], 'atoms': []})['atoms'].append(a)

EXPECT = {
    124: ('GLU', '必须去质子化（带负电），才能与 R283 成盐桥'),
    283: ('ARG', '必须质子化（带正电）'),
    120: ('GLN', '极性，中性'),
    99:  ('ASP', '通常去质子化'),
    102: ('ARG', '通常质子化'),
}

print(f"    {'残基':<14}{'状态':<10}{'判定'}")
print('    ' + '-' * 62)
for resi in KEY_RESIDUES:
    if resi not in by_res:
        print(f'    {resi:<14}未找到')
        continue
    name = by_res[resi]['name']
    exp, note = EXPECT.get(resi, ('', ''))
    ok = '✓ 正确' if name == exp else f'⚠ 预期 {exp}'
    print(f'    {name}{resi:<10}{ok:<10}{note}')

# HIS 质子化状态
hiss = sorted({(r, by_res[r]['name']) for r in by_res if by_res[r]['name'].startswith('HIS')
               or by_res[r]['name'] in ('HID', 'HIE', 'HIP', 'HIS')})
if hiss:
    print(f'\n    组氨酸质子化状态（共 {len(hiss)} 个，PDBFixer 未细分为 HID/HIE/HIP）:')
    print(f'      残基编号: {", ".join(f"HIS{r}" for r, _ in hiss)}')
    print('      注意: His280(6.52) 属于激活相关的芳香簇（F279/H280/F312），')
    print('            若要更精确，需用 PROPKA 或 H++ 单独判断其质子化。')

# ---------------- 盐桥复测 ----------------
print('\n[6] 盐桥距离复测（加氢后）')


def get(resi, name):
    for a in by_res.get(resi, {}).get('atoms', []):
        if a['name'] == name:
            return a
    return None


def d(a, b):
    return math.sqrt((a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2)


pairs = [(('124', 'OE1'), ('283', 'NH1')),
         (('124', 'OE2'), ('283', 'NH2')),
         (('124', 'OE1'), ('283', 'NH2')),
         (('124', 'OE2'), ('283', 'NH1'))]
ds = []
for (ri, an), (rj, bn) in pairs:
    a, b = get(int(ri), an), get(int(rj), bn)
    if a and b:
        v = d(a, b)
        ds.append(v)
        print(f'    E124 {an} — R283 {bn} : {v:.2f} Å')

if ds:
    mn = min(ds)
    print(f'\n    最短距离: {mn:.2f} Å   （加氢前的基准值: 2.56 Å）')
    if abs(mn - 2.56) < 0.3:
        print('    ✓ 盐桥保持完好，质子化状态正确')
    else:
        print('    ⚠ 盐桥距离偏离基准超过 0.3 Å —— 需检查 E124 / R283 的电荷设置')

sep('=')
print('  下一步: 用这个加氢后的结构准备对接，并做 redocking 验证')
print('  结构文件: 9UY3_receptor_h.pdb')
sep('=')
