# -*- coding: utf-8 -*-
"""XGBoost / LightGBM vs 当前 GradientBoosting 对比"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import spearmanr

try:
    from xgboost import XGBRegressor
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
try:
    from lightgbm import LGBMRegressor
    HAS_LGB = True
except ImportError:
    HAS_LGB = False

RDLogger.DisableLog('rdApp.*')
SEEDS = [0, 1, 2, 3, 4]

print('=== 核心库版本（确认未被降级）===')
import rdkit, sklearn
print(f'  rdkit {rdkit.__version__} / numpy {np.__version__} / pandas {pd.__version__} / sklearn {sklearn.__version__}')
print(f'  xgboost {"已装" if HAS_XGB else "未装"}   lightgbm {"已装" if HAS_LGB else "未装"}')
print()

# ---------- 数据 ----------
df = pd.read_csv('ghsr1a_binding_clean.csv')
df = df[~df['peptide_like']].reset_index(drop=True)
df['mol'] = df['canonical_smiles'].apply(Chem.MolFromSmiles)
df = df[df['mol'].notna()].reset_index(drop=True)
n = len(df); y = df['pAffinity'].values; mols = df['mol'].tolist()
print(f'数据集: {n} 个 GHSR1a 小分子')

DESC = {
    'MW': Descriptors.MolWt, 'logP': Crippen.MolLogP,
    'TPSA': rdMolDescriptors.CalcTPSA, 'HBD': rdMolDescriptors.CalcNumHBD,
    'HBA': rdMolDescriptors.CalcNumHBA, 'RotB': rdMolDescriptors.CalcNumRotatableBonds,
    'AromR': rdMolDescriptors.CalcNumAromaticRings, 'Rings': rdMolDescriptors.CalcNumRings,
    'SatR': rdMolDescriptors.CalcNumSaturatedRings,
    'HeavyAtoms': Descriptors.HeavyAtomCount,
    'FracCsp3': rdMolDescriptors.CalcFractionCSP3, 'QED': QED.qed,
    'MolMR': Crippen.MolMR, 'BertzCT': Descriptors.BertzCT,
    'FpDensity': Descriptors.FpDensityMorgan1,
}
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.array([gen.GetFingerprint(m).ToList() for m in mols], dtype=np.float32)
desc = np.nan_to_num(np.array([[f(m) for f in DESC.values()] for m in mols], dtype=np.float32))
X = np.hstack([morgan, desc])
print(f'特征维度: {X.shape[1]}（Morgan2048 + 15 描述符）\n')

# ---------- 切分 ----------
scaf = defaultdict(list)
for i, s in enumerate(df['scaffold'].values):
    scaf[s].append(i)


def scaffold_split(seed=0):
    rng = np.random.RandomState(seed)
    gs = list(scaf.values())
    rng.shuffle(gs)
    gs.sort(key=len, reverse=True)
    nt = int(n * 0.8)
    tr, te = [], []
    for g in gs:
        (tr if len(tr) + len(g) <= nt else te).extend(g)
    return np.array(tr), np.array(te)


# ---------- 模型 ----------
MODELS = {
    'GradientBoosting (当前最佳)':
        lambda s: GradientBoostingRegressor(n_estimators=300, random_state=s),
    'RandomForest (对照)':
        lambda s: RandomForestRegressor(n_estimators=500, random_state=s, n_jobs=-1),
}

if HAS_XGB:
    MODELS['XGBoost (默认)'] = lambda s: XGBRegressor(
        n_estimators=300, random_state=s, n_jobs=-1, verbosity=0)
    MODELS['XGBoost (调参)'] = lambda s: XGBRegressor(
        n_estimators=600, learning_rate=0.03, max_depth=5,
        subsample=0.8, colsample_bytree=0.7,
        reg_alpha=0.1, reg_lambda=1.0,
        min_child_weight=2, random_state=s, n_jobs=-1, verbosity=0)

if HAS_LGB:
    MODELS['LightGBM (默认)'] = lambda s: LGBMRegressor(
        n_estimators=300, random_state=s, n_jobs=-1, verbose=-1)
    MODELS['LightGBM (调参)'] = lambda s: LGBMRegressor(
        n_estimators=600, learning_rate=0.03, num_leaves=31,
        subsample=0.8, colsample_bytree=0.7,
        reg_alpha=0.1, reg_lambda=1.0,
        min_child_samples=5, random_state=s, n_jobs=-1, verbose=-1)


def evaluate(mfn):
    r2s, maes, rmses, rhos = [], [], [], []
    for s in SEEDS:
        tr, te = scaffold_split(s)
        m = mfn(s)
        m.fit(X[tr], y[tr])
        p = m.predict(X[te])
        r2s.append(r2_score(y[te], p))
        maes.append(mean_absolute_error(y[te], p))
        rmses.append(np.sqrt(mean_squared_error(y[te], p)))
        rhos.append(spearmanr(y[te], p).correlation)
    return np.mean(r2s), np.std(r2s), np.mean(maes), np.mean(rmses), np.mean(rhos)


print('=' * 82)
print('  Scaffold Split · 5 seeds · 模型对比')
print('=' * 82)
print(f"  {'模型':<28}{'R²':>15}{'MAE':>8}{'RMSE':>8}{'Spearman':>10}")
print('  ' + '-' * 78)

rows = []
for name, mfn in MODELS.items():
    r2, r2s, mae, rmse, rho = evaluate(mfn)
    rows.append((name, r2, r2s, mae, rmse, rho))
    print(f'  {name:<26}{r2:>10.3f}±{r2s:<4.2f}{mae:>8.3f}{rmse:>8.3f}{rho:>10.3f}')

print('  ' + '-' * 78)

R = pd.DataFrame(rows, columns=['模型', 'R2', 'R2波动', 'MAE', 'RMSE', 'Spearman'])
R = R.sort_values('R2', ascending=False).reset_index(drop=True)
best = R.iloc[0]
base = R[R['模型'].str.contains('GradientBoosting')].iloc[0]

print(f'\n  最佳: {best["模型"]}')
print(f'    R² = {best["R2"]:.3f}   MAE = {best["MAE"]:.3f}   Spearman = {best["Spearman"]:.3f}')
print(f'\n  相对当前 GradientBoosting 的提升:')
print(f'    R²       {base["R2"]:.3f} → {best["R2"]:.3f}   ({(best["R2"]-base["R2"])/base["R2"]*100:+.1f}%)')
print(f'    MAE      {base["MAE"]:.3f} → {best["MAE"]:.3f}   ({(base["MAE"]-best["MAE"])/base["MAE"]*100:+.1f}% 误差下降)')
print(f'    Spearman {base["Spearman"]:.3f} → {best["Spearman"]:.3f}')
print(f'\n  基线（永远预测均值）MAE = {np.abs(y - y.mean()).mean():.3f}')
print(f'  最佳模型相对基线提升 {(np.abs(y-y.mean()).mean()-best["MAE"])/np.abs(y-y.mean()).mean()*100:.0f}%')

R.to_csv('boosting_comparison.csv', index=False)
print('\n  已保存: boosting_comparison.csv')
