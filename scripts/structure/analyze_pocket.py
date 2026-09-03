# -*- coding: utf-8 -*-
"""
GHSR1a 结合口袋分析
算出 anamorelin (A1ESD) 被哪些残基包围，以及各残基属于哪根跨膜螺旋。
纯标准库实现，不依赖 mdtraj。
"""
import os, math, requests
from collections import defaultdict

D = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(D, '9UY3_receptor.pdb')
LIG = os.path.join(D, '9UY3_ligand_A1ESD.pdb')

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


def read_pdb(path):
    atoms = []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            atoms.append({
                'name': line[12:16].strip(),
                'resn': line[17:20].strip(),
                'chain': line[21],
                'resi': int(line[22:26]),
                'x': float(line[30:38]), 'y': float(line[38:46]),
                'z': float(line[46:54]),
                'elem': line[76:78].strip() or line[12:16].strip()[0],
            })
    return atoms


def dist(a, b):
    return math.sqrt((a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2)


print('=' * 74)
print('  GHSR1a 结合口袋分析 · 9UY3 + anamorelin')
print('=' * 74)

# ---------- 1. 从 UniProt 取跨膜区边界 ----------
print('\n[1] 获取跨膜螺旋边界 (UniProt Q92847)')
r = requests.get('https://rest.uniprot.org/uniprotkb/Q92847.json', timeout=60).json()
tms = []
for f in r.get('features', []):
    if f['type'] == 'Transmembrane':
        s = f['location']['start']['value']
        e = f['location']['end']['value']
        name = f.get('description', '').split('Name=')[-1]
        tms.append((s, e, name))
tms.sort()
for s, e, n in tms:
    print(f'    TM{n}  {s:>3} – {e:<3}')


def which_tm(resi):
    for s, e, n in tms:
        if s <= resi <= e:
            return f'TM{n}'
    # 环区。GPCR 拓扑: N 端在胞外, 依次
    #   TM1 →ICL1→ TM2 →ECL1→ TM3 →ICL2→ TM4 →ECL2→ TM5 →ICL3→ TM6 →ECL3→ TM7
    # 即 i 为偶数(TM1-TM2, TM3-TM4, TM5-TM6)是胞内环, 奇数(TM2-TM3, TM4-TM5, TM6-TM7)是胞外环
    for i in range(len(tms) - 1):
        s, e, _ = tms[i]
        nxt = tms[i + 1][0]
        if e < resi < nxt:
            if i % 2 == 0:
                return f'ICL{i // 2 + 1}'
            else:
                return f'ECL{(i + 1) // 2}'
    if resi < tms[0][0]:
        return 'N-term'
    return 'C-term'


# ---------- 2. 读结构 ----------
rec = read_pdb(REC)
lig = read_pdb(LIG)
print(f'\n[2] 读取结构')
print(f'    受体 {len(rec)} 原子 / {len({a["resi"] for a in rec})} 残基')
print(f'    配体 {len(lig)} 原子 (anamorelin)')

# ---------- 3. 算最小距离 ----------
print(f'\n[3] 计算每个残基到配体的最小距离')
res_min = {}
res_atom = {}
for a in rec:
    if a['elem'] == 'H':
        continue
    d = min(dist(a, b) for b in lig)
    key = a['resi']
    if key not in res_min or d < res_min[key]:
        res_min[key] = d
        res_atom[key] = a['name']

# ---------- 4. 输出口袋残基 ----------
CUT = 5.0
pocket = sorted([(k, v) for k, v in res_min.items() if v <= CUT], key=lambda x: x[1])
res_name = {a['resi']: a['resn'] for a in rec}

print(f'\n[4] 配体 {CUT} Å 内的残基（共 {len(pocket)} 个）')
print(f"    {'残基':<12}{'TM 归属':<10}{'最小距离':>9}   接触原子")
print('    ' + '-' * 52)
for resi, d in pocket:
    rn = res_name[resi]
    one = AA3TO1.get(rn, 'X')
    print(f'    {rn}{resi:<8}{which_tm(resi):<10}{d:>7.2f} Å   {res_atom[resi]}')

# ---------- 5. 按螺旋汇总 ----------
print(f'\n[5] 按结构区域汇总')
by_region = defaultdict(list)
for resi, d in pocket:
    by_region[which_tm(resi)].append(f'{AA3TO1.get(res_name[resi],"X")}{resi}')
print(f"    {'区域':<10}{'残基数':>7}   残基")
print('    ' + '-' * 62)
order = [f'TM{n}' for n in range(1, 8)] + ['ECL1', 'ECL2', 'ECL3', 'ICL1', 'ICL2', 'ICL3']
for reg in order:
    if reg in by_region:
        lst = by_region[reg]
        print(f'    {reg:<10}{len(lst):>7}   {", ".join(lst)}')
for reg in sorted(by_region):
    if reg not in order:
        lst = by_region[reg]
        print(f'    {reg:<10}{len(lst):>7}   {", ".join(lst)}')

# ---------- 6. 关键提示 ----------
print(f'\n[6] 读这张表的要点')
charged = [(resi, res_name[resi]) for resi, d in pocket
           if res_name[resi] in ('ASP', 'GLU', 'ARG', 'LYS', 'HIS')]
if charged:
    print('    带电荷残基（最影响结合，也最容易因质子化状态出错）:')
    for resi, rn in charged:
        print(f'      {rn}{resi}  ({which_tm(resi)})  距离 {res_min[resi]:.2f} Å')
arom = [(resi, res_name[resi]) for resi, d in pocket
        if res_name[resi] in ('PHE', 'TYR', 'TRP', 'HIS')]
if arom:
    print('    芳香残基（π 堆积 / 阳离子-π 作用）:')
    for resi, rn in arom:
        print(f'      {rn}{resi}  ({which_tm(resi)})  距离 {res_min[resi]:.2f} Å')

print(f'\n    保存为 pocket_residues.txt 备查')
with open(os.path.join(D, 'pocket_residues.txt'), 'w') as f:
    f.write('GHSR1a pocket residues within 5 A of anamorelin (9UY3)\n\n')
    for resi, d in pocket:
        f.write(f'{res_name[resi]}{resi}\t{which_tm(resi)}\t{d:.2f}\n')
