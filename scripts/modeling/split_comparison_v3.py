# -*- coding: utf-8 -*-
"""切分严格度连续光谱 + 大样本分层诊断"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from collections import defaultdict
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')

SEEDS = [0, 1, 2, 3, 4]

df = pd.read_csv('ghsr1a_binding_clean.csv')
df = df[~df['peptide_like']].reset_index(drop=True)
df['mol'] = df['canonical_smiles'].apply(Chem.MolFromSmiles)
df = df[df['mol'].notna()].reset_index(drop=True)
n = len(df)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps_list = [gen.GetFingerprint(m) for m in df['mol']]
X = np.array([f.ToList() for f in fps_list], dtype=np.float32)
y = df['pAffinity'].values
print(f'数据集: GHSR1a 小分子子集  n = {n}\n')

# 预计算距离矩阵（Butina 各 cutoff 共用）
print('计算相似度矩阵 ...')
dists = []
for i in range(n):
    s = DataStructs.BulkTanimotoSimilarity(fps_list[i], fps_list[:i])
    dists.extend([1 - x for x in s])
print(f'  {n*(n-1)//2} 个分子对完成\n')


def cluster_at(cutoff):
    return Butina.ClusterData(dists, n, cutoff, isDistData=True)


def group_split(groups, frac=0.8, seed=0):
    rng = np.random.RandomState(seed)
    gs = list(groups)
    rng.shuffle(gs)
    gs.sort(key=len, reverse=True)
    n_train = int(n * frac)
    tr, te = [], []
    for g in gs:
        if len(tr) + len(g) <= n_train:
            tr.extend(g)
        else:
            te.extend(g)
    return np.array(tr), np.array(te)


def random_split(frac=0.8, seed=0):
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    k = int(n * frac)
    return idx[:k], idx[k:]


def max_sim_to_train(te, tr, X):
    T, R = X[te], X[tr]
    out = []
    for i in range(0, len(T), 128):
        b = T[i:i + 128]
        inter = b @ R.T
        union = b.sum(1, keepdims=True) + R.sum(1, keepdims=True).T - inter
        out.append(np.where(union > 0, inter / union, 0).max(axis=1))
    return np.concatenate(out)


def run(tr, te, seed):
    m = RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=-1)
    m.fit(X[tr], y[tr])
    p = m.predict(X[te])
    return dict(r2=r2_score(y[te], p),
                mae=mean_absolute_error(y[te], p),
                rmse=np.sqrt(mean_squared_error(y[te], p)),
                rho=spearmanr(y[te], p).correlation,
                sim=np.median(max_sim_to_train(te, tr, X)))


def eval_split(getter, label, scenario):
    rows = [run(*getter(s), seed=s) for s in SEEDS]
    r = pd.DataFrame(rows)
    return dict(label=label, scenario=scenario, r2=r['r2'].mean(), r2sd=r['r2'].std(),
                mae=r['mae'].mean(), rho=r['rho'].mean(), sim=r['sim'].mean())


# 骨架分组
scaf_groups = defaultdict(list)
for i, s in enumerate(df['scaffold'].values):
    scaf_groups[s].append(i)

results = []
results.append(eval_split(
    lambda s: random_split(seed=s), '随机切分', '同系列补测数据（过度乐观）'))
results.append(eval_split(
    lambda s: group_split(scaf_groups.values(), seed=s), 'Scaffold Split',
    '新骨架、同化学类型'))

for cutoff in [0.35, 0.45, 0.55, 0.65]:
    cl = cluster_at(cutoff)
    results.append(eval_split(
        lambda s, cl=cl: group_split(cl, seed=s), f'Butina {cutoff:.2f}',
        '新化学系列' if cutoff <= 0.45 else '跨化学类型（不现实）'))

R = pd.DataFrame(results)

print('=' * 92)
print('  切分严格度光谱')
print('=' * 92)
print(f"  {'切分方式':<16}{'最大相似度':>10}{'R²':>14}{'MAE':>8}{'Spearman':>10}  {'对应真实场景'}")
print('  ' + '-' * 88)
for _, r in R.iterrows():
    print(f"  {r['label']:<14}{r['sim']:>10.3f}{r['r2']:>10.3f}±{r['r2sd']:.2f}"
          f"{r['mae']:>8.3f}{r['rho']:>10.3f}  {r['scenario']}")

base = np.abs(y - y.mean()).mean()
print(f"\n  参考基线（永远预测训练集均值）MAE = {base:.3f}")
print('  高于此值即模型无效\n')

# ---------------- 大样本分层诊断 ----------------
print('=' * 92)
print('  分层诊断 · 模型在多"新"的分子上才开始失效（10 次重复）')
print('=' * 92)
rows = []
for s in range(10):
    idx = np.random.RandomState(100 + s).permutation(n)
    k = int(n * 0.8)
    tr, te = idx[:k], idx[k:]
    m = RandomForestRegressor(n_estimators=300, random_state=s, n_jobs=-1)
    m.fit(X[tr], y[tr])
    p = m.predict(X[te])
    rows.append(pd.DataFrame({'sim': max_sim_to_train(te, tr, X),
                              'err': np.abs(y[te] - p)}))
D = pd.concat(rows, ignore_index=True)
print(f"  累计 {len(D)} 个测试分子\n")
print(f"  {'相似度区间':<18}{'分子数':>8}{'占比':>8}{'MAE':>8}   判读")
print('  ' + '-' * 62)
for lo, hi, note in [(0.0, 0.5, '基本无关'), (0.5, 0.6, '远亲'), (0.6, 0.7, '同一化学系列'),
                     (0.7, 0.85, '近亲'), (0.85, 1.01, '高度相似')]:
    sub = D[(D['sim'] >= lo) & (D['sim'] < hi)]
    if len(sub) == 0:
        continue
    mae = sub['err'].mean()
    verdict = '模型失效' if mae > base else ('勉强' if mae > base * 0.85 else '可靠')
    print(f"  {lo:.2f} – {min(hi,1.0):.2f}{'':<8}{len(sub):>8}{len(sub)/len(D)*100:>7.1f}%"
          f"{mae:>8.3f}   {note} · {verdict}")

print(f"\n  结论：模型失效的临界点在与训练集相似度约 "
      f"{D[D['err'] > base]['sim'].median():.2f} 处")

R.to_csv('split_comparison_results.csv', index=False)
print('\n  结果已保存: split_comparison_results.csv')
