# -*- coding: utf-8 -*-
"""s1_18d_L2_pot.py — L2 起源裁决 (POT 纯 numpy 版, 第7次也是最后一次)
用 POT 的 Sinkhorn 直接算 OT 传输矩阵, 不依赖 jax/moscot/cellrank。
链式: 24h --OT--> 7d --OT--> 6mo
富集比 R = mean(行和到 6mo p16+ | 急性 p21+) / mean(... | p21−)
"""
import os, time
import numpy as np, pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy import sparse
import ot

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
INPD = os.path.join(SEN, "results/s1_L2_origin")
SEED = 20260925
np.random.seed(SEED)
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

DAY = {"Ctrl": 0.0, "24h": 1.0, "7d": 7.0, "6mo": 180.0}
CAP = 4000  # POT Sinkhorn 对 4k×4k 足够快

t0 = time.time()
log("[L2d] 读入 ...")
X = mmread(os.path.join(INPD, "counts.mtx")).T.tocsr()
genes = [l.strip() for l in open(os.path.join(INPD, "genes.txt"))]
cells = pd.read_csv(os.path.join(INPD, "cells.csv"), index_col="cell_id")
adata = ad.AnnData(X=X.astype(np.float32), obs=cells, var=pd.DataFrame(index=pd.Index(genes)))
adata.obs["day"] = adata.obs["time_point"].map(DAY).astype(float)

# 截断
keep = []
for tp_name, idx in adata.obs.groupby("time_point").groups.items():
    idx = np.asarray(idx)
    if len(idx) > CAP: idx = np.random.choice(idx, CAP, replace=False)
    keep += list(idx)
adata = adata[sorted(keep)].copy()
log(f"[L2d] 截断后: {dict(adata.obs.time_point.value_counts())}")

# 预处理 → 潜入空间
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, svd_solver="arpack", random_state=SEED)
try:
    sc.external.pp.harmony_integrate(adata, key="lib", basis="X_pca",
                                     adjusted_basis="X_pca_h", max_iter_harmony=20)
    emb = adata.obsm["X_pca_h"]
except:
    emb = adata.obsm["X_pca"]
log(f"[L2d] 潜入空间 {emb.shape} ({time.time()-t0:.0f}s)")

# 状态旗标
acute = (adata.obs["day"] == 1.0).values
term_p16 = ((adata.obs["day"] == 180.0) & (adata.obs["Cdkn2a"] > 0)).values
p21pos_ac = (adata.obs["Cdkn1a"] > 0).values[acute]
log(f"[L2d] 24h n={acute.sum()} p21+={p21pos_ac.sum()} | 6mo p16+ n={term_p16.sum()}")

# ===== OT 传输 =====
def sinkhorn_ot(X_src, X_tgt, reg=0.005):
    """POT Sinkhorn: 返回归一化传输矩阵 (n_src × n_tgt), 行和=1"""
    n, m = len(X_src), len(X_tgt)
    # 代价 = 欧氏距离
    M = ot.dist(X_src, X_tgt, metric="euclidean")
    M /= (M.max() + 1e-12)
    # 均匀分布 (不做生长/死亡建模 — 只看状态流)
    a_w = np.full(n, 1/n); b_w = np.full(m, 1/m)
    T = ot.sinkhorn(a_w, b_w, M, reg=reg, numItermax=5000, stopThr=1e-9)
    return T

idx_24h = np.where(acute)[0]
idx_7d  = np.where(adata.obs["day"] == 7.0)[0]
idx_6mo = np.where(adata.obs["day"] == 180.0)[0]

log(f"\n[L2d] OT 24h({len(idx_24h)}) → 7d({len(idx_7d)}) ...")
P1 = sinkhorn_ot(emb[idx_24h], emb[idx_7d])
log(f"  完成 ({time.time()-t0:.0f}s)")

log(f"[L2d] OT 7d({len(idx_7d)}) → 6mo({len(idx_6mo)}) ...")
P2 = sinkhorn_ot(emb[idx_7d], emb[idx_6mo])
log(f"  完成 ({time.time()-t0:.0f}s)")

# 链式: 24h → (7d) → 6mo
Pch = P1 @ P2  # n_24h × n_6mo
# 提取: 急性 p21+ vs p21− 传输到 6mo p16+ 的质量
idx_6mo_p16_in6mo = term_p16[idx_6mo]  # 在 idx_6mo 子集内的 p16+ 位置
mass = Pch[:, idx_6mo_p16_in6mo].sum(axis=1)  # 每 24h 细胞传到 6mo p16+ 的总质量

R = mass[p21pos_ac].mean() / max(mass[~p21pos_ac].mean(), 1e-12)
m_pos, m_neg = mass[p21pos_ac].mean(), mass[~p21pos_ac].mean()
log(f"\n[L2d] p21+→p16+ 质量: {m_pos:.5f} vs p21−: {m_neg:.5f}")
log(f"[L2d] 富集比 R = {R:.3f}")

# 置换 (1000 次)
rng = np.random.default_rng(SEED); cnt = 0
for _ in range(1000):
    perm = rng.permutation(p21pos_ac)
    r = mass[perm].mean() / max(mass[~perm].mean(), 1e-12)
    if r >= R: cnt += 1
p_perm = (cnt + 1) / 1001
log(f"[L2d] 置换 p = {p_perm:.4f} (n=1000)")

# 副法: 单步 24h→7d
idx_7d_p16 = ((adata.obs["day"] == 7.0) & (adata.obs["Cdkn2a"] > 0)).values[idx_7d]
mass_7d = P1[:, idx_7d_p16].sum(axis=1)
R7 = mass_7d[p21pos_ac].mean() / max(mass_7d[~p21pos_ac].mean(), 1e-12)
log(f"[L2d] 24h→7d p16+ 单步 R = {R7:.3f} (7d p16+ n={idx_7d_p16.sum()})")

# 判定
log("\n" + "=" * 72)
verdict = "延续" if (R >= 2 and p_perm < 0.05) else ("独立起源" if 0.5 <= R <= 2 else "中间/不可判")
log(f"判定: 链式 R={R:.3f}, p_perm={p_perm:.4f} → {verdict}")
log(f"       单步 R={R7:.3f}")

log(f"\n[L2d] DONE ({time.time()-t0:.0f}s)")
with open(os.path.join(INPD, "L2_results.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
pd.DataFrame([{"method": "POT_sinkhorn", "R_chain": R, "p_perm": p_perm, "R_step7d": R7,
               "m_pos": m_pos, "m_neg": m_neg}]).to_csv(os.path.join(INPD, "L2_enrichment.csv"), index=False)
