# -*- coding: utf-8 -*-
"""QSAR 基线实验：特征 × 模型 网格对比（scaffold split）"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys
from sklearn.ensemble import (RandomForestRegressor, ExtraTreesRegressor,
                              GradientBoostingRegressor)
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')

SEEDS = [0, 1, 2, 3, 4]

# ---------- 数据 ----------
df = pd.read_csv('ghsr1a_binding_clean.csv')
df = df[~df['peptide_like']].reset_index(drop=True)
df['mol'] = df['canonical_smiles'].apply(Chem.MolFromSmiles)
df = df[df['mol'].notna()].reset_index(drop=True)
n = len(df)
y = df['pAffinity'].values
print(f'数据集: {n} 个 GHSR1a 小分子')
print(f'pAffinity: 均值 {y.mean():.2f}  std {y.std():.2f}\n')

mols = df['mol'].tolist()

# ---------- 特征 ----------
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
gen1024 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

DESC_FUNCS = {
    'MW': Descriptors.MolWt,
    'logP': Crippen.MolLogP,
    'TPSA': rdMolDescriptors.CalcTPSA,
    'HBD': rdMolDescriptors.CalcNumHBD,
    'HBA': rdMolDescriptors.CalcNumHBA,
    'RotB': rdMolDescriptors.CalcNumRotatableBonds,
    'AromR': rdMolDescriptors.CalcNumAromaticRings,
    'Rings': rdMolDescriptors.CalcNumRings,
    'SatR': rdMolDescriptors.CalcNumSaturatedRings,
    'HeavyAtoms': Descriptors.HeavyAtomCount,
    'FracCsp3': rdMolDescriptors.CalcFractionCSP3,
    'QED': QED.qed,
    'MolMR': Crippen.MolMR,
    'BertzCT': Descriptors.BertzCT,
    'FpDensity': Descriptors.FpDensityMorgan1,
}

print('计算特征 ...')
morgan2048 = np.array([gen.GetFingerprint(m).ToList() for m in mols], dtype=np.float32)
morgan1024 = np.array([gen1024.GetFingerprint(m).ToList() for m in mols], dtype=np.float32)
maccs = np.array([list(MACCSkeys.GenMACCSKeys(m))[1:] for m in mols], dtype=np.float32)
desc = np.array([[f(m) for f in DESC_FUNCS.values()] for m in mols], dtype=np.float32)
desc = np.nan_to_num(desc, nan=0.0, posinf=0.0, neginf=0.0)
print(f'  Morgan2048 {morgan2048.shape}  MACCS {maccs.shape}  描述符 {desc.shape}\n')

FEATURES = {
    'Morgan2048': morgan2048,
    'Morgan1024': morgan1024,
    'MACCS(167)': maccs,
    '描述符(15)': desc,
    'Morgan+描述符': np.hstack([morgan2048, desc]),
    'Morgan+MACCS+描述符': np.hstack([morgan2048, maccs, desc]),
}

MODELS = {
    'RandomForest': lambda s: RandomForestRegressor(n_estimators=500, random_state=s, n_jobs=-1),
    'ExtraTrees':   lambda s: ExtraTreesRegressor(n_estimators=500, random_state=s, n_jobs=-1),
    'GradBoost':    lambda s: GradientBoostingRegressor(n_estimators=300, random_state=s),
}

# ---------- scaffold split ----------
scaf_groups = defaultdict(list)
for i, s in enumerate(df['scaffold'].values):
    scaf_groups[s].append(i)


def scaffold_split(seed):
    rng = np.random.RandomState(seed)
    gs = list(scaf_groups.values())
    rng.shuffle(gs)
    gs.sort(key=len, reverse=True)
    n_train = int(n * 0.8)
    tr, te = [], []
    for g in gs:
        if len(tr) + len(g) <= n_train:
            tr.extend(g)
        else:
            te.extend(g)
    return np.array(tr), np.array(te)


def random_split(seed):
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    k = int(n * 0.8)
    return idx[:k], idx[k:]


# ---------- 实验 ----------
def evaluate(X, model_fn, splitter):
    r2s, maes, rmses, rhos = [], [], [], []
    for s in SEEDS:
        tr, te = splitter(s)
        m = model_fn(s)
        m.fit(X[tr], y[tr])
        p = m.predict(X[te])
        r2s.append(r2_score(y[te], p))
        maes.append(mean_absolute_error(y[te], p))
        rmses.append(np.sqrt(mean_squared_error(y[te], p)))
        rhos.append(spearmanr(y[te], p).correlation)
    return (np.mean(r2s), np.std(r2s), np.mean(maes),
            np.mean(rmses), np.mean(rhos))


print('=' * 84)
print('  Scaffold Split 下的特征 × 模型对比')
print('=' * 84)
print(f"  {'特征':<22}{'模型':<14}{'R²':>13}{'MAE':>8}{'RMSE':>8}{'Spearman':>10}")
print('  ' + '-' * 80)

rows = []
for fname, X in FEATURES.items():
    for mname, mfn in MODELS.items():
        r2, r2s, mae, rmse, rho = evaluate(X, mfn, scaffold_split)
        rows.append((fname, mname, r2, mae, rmse, rho))
        print(f'  {fname:<20}{mname:<14}{r2:>8.3f}±{r2s:<4.2f}{mae:>8.3f}{rmse:>8.3f}{rho:>10.3f}')

R = pd.DataFrame(rows, columns=['特征', '模型', 'R2', 'MAE', 'RMSE', 'Spearman'])
best = R.loc[R['R2'].idxmax()]
print('  ' + '-' * 80)
print(f"\n  最佳组合: {best['特征']} + {best['模型']}")
print(f"    R² = {best['R2']:.3f}   MAE = {best['MAE']:.3f}   Spearman = {best['Spearman']:.3f}")

# ---------- 同一最佳组合下：随机切分 vs scaffold ----------
print('\n' + '=' * 84)
print('  最佳组合的切分方式对比')
print('=' * 84)
Xb = FEATURES[best['特征']]
mfn = MODELS[best['模型']]
for sname, splitter in [('随机切分', random_split), ('Scaffold Split', scaffold_split)]:
    r2, r2s, mae, rmse, rho = evaluate(Xb, mfn, splitter)
    print(f'  {sname:<16} R² = {r2:.3f}±{r2s:.2f}   MAE = {mae:.3f}   Spearman = {rho:.3f}')

# ---------- 基线 ----------
print(f'\n  参考基线（永远预测训练集均值）MAE = {np.abs(y - y.mean()).mean():.3f}')

R.to_csv('qsar_feature_model_grid.csv', index=False)
print('\n  结果已保存: qsar_feature_model_grid.csv')
