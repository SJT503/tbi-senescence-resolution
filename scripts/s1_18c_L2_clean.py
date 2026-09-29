# -*- coding: utf-8 -*-
"""s1_18c_L2_clean.py — L2 起源裁决 (干净重写版, 第4次)
主法: moscot TemporalProblem — 链式 24h→7d→6mo
富集比 R = mean(行和到 6mo p16+ | 急性 p21+) / mean(... | 急性 p21−)
副法: 单步 24h→7d
判据(预登记 L2_PROTOCOL.md): R>=2 & p_perm<0.05 → 延续; R in 0.5-2 & ns → 独立起源
"""
import os, time
import numpy as np, pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy import sparse

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
INPD = os.path.join(SEN, "results/s1_L2_origin")
SEED = 20260925
np.random.seed(SEED)
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

DAY = {"Ctrl": 0.0, "24h": 1.0, "7d": 7.0, "6mo": 180.0}
CAP_PER_TP = 6000

t0 = time.time()
log("[L2c] 读入 ...")
X = mmread(os.path.join(INPD, "counts.mtx")).T.tocsr()
genes = [l.strip() for l in open(os.path.join(INPD, "genes.txt"))]
cells = pd.read_csv(os.path.join(INPD, "cells.csv"), index_col="cell_id")
assert X.shape == (len(cells), len(genes))
adata = ad.AnnData(X=X.astype(np.float32), obs=cells, var=pd.DataFrame(index=pd.Index(genes)))
adata.obs["day"] = adata.obs["time_point"].map(DAY).astype(float)
log(f"[L2c] {adata.n_obs} x {adata.n_vars}")

# 截断
keep = []
for tp_name, idx in adata.obs.groupby("time_point").groups.items():
    idx = np.asarray(idx)
    if len(idx) > CAP_PER_TP:
        idx = np.random.choice(idx, CAP_PER_TP, replace=False)
        log(f"  截断 {tp_name}→{CAP_PER_TP}")
    keep += list(idx)
adata = adata[sorted(keep)].copy()

# 预处理
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, svd_solver="arpack", random_state=SEED)
try:
    sc.external.pp.harmony_integrate(adata, key="lib", basis="X_pca",
                                     adjusted_basis="X_pca_harmony", max_iter_harmony=20)
    log("[L2c] Harmony OK")
except Exception as e:
    log(f"[L2c] Harmony skip: {e}")
sc.pp.neighbors(adata, use_rep="X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca",
                n_neighbors=15, random_state=SEED)
log(f"[L2c] 预处理完成 ({time.time()-t0:.0f}s)")

# 状态旗标
acute = (adata.obs["day"] == 1.0).values
term_p16 = ((adata.obs["day"] == 180.0) & (adata.obs["Cdkn2a"] > 0)).values
a = (adata.obs["Cdkn1a"] > 0).values[acute]
log(f"[L2c] 24h n={acute.sum()} p21+={a.sum()} | 6mo p16+ n={term_p16.sum()}")

results = {}

# ===== 主法: moscot 链式 =====
log("\n[L2c] === moscot TemporalProblem ===")
try:
    from moscot.problems import TemporalProblem
    tp = TemporalProblem(adata)
    tp = tp.prepare(time_key="day")
    tp = tp.solve(epsilon=0.005, tau_a=0.99, tau_b=0.99)  # ⚠️ 无 seed 参数
    log(f"[L2c] moscot solve 完成 ({time.time()-t0:.0f}s)")

    # 链式 24h→6mo (⚠️ moscot 返回 jax array → 需转 scipy sparse)
    P1 = sparse.csr_matrix(np.asarray(tp[1.0, 7.0].solution.transport_matrix))    # 24h x 7d
    P2 = sparse.csr_matrix(np.asarray(tp[7.0, 180.0].solution.transport_matrix))  # 7d x 6mo
    P1.eliminate_zeros(); P2.eliminate_zeros()
    Pch = (P1 @ P2).tocsr()
    rs = np.asarray(Pch.sum(axis=1)).ravel(); rs[rs == 0] = 1
    Pch = sparse.diags(1/rs) @ Pch
    idx_ac = np.where(acute)[0]; idx_tm = np.where(term_p16)[0]
    mass = np.asarray(Pch[np.ix_(idx_ac, idx_tm)].sum(axis=1)).ravel()
    R = mass[a].mean() / max(mass[~a].mean(), 1e-12)
    # 置换
    rng = np.random.default_rng(SEED); cnt = 0
    for _ in range(1000):
        perm = rng.permutation(a)
        r = mass[perm].mean() / max(mass[~perm].mean(), 1e-12)
        if r >= R: cnt += 1
    p_perm = (cnt + 1) / 1001
    results["moscot_chain"] = dict(R=R, p_perm=p_perm, m_pos=mass[a].mean(), m_neg=mass[~a].mean())
    log(f"  p21+→6mo p16+ 质量: {mass[a].mean():.4f} vs p21−: {mass[~a].mean():.4f}")
    log(f"  链式富集比 R = {R:.3f} | 置换 p = {p_perm:.4f}")
    np.save(os.path.join(INPD, "_L2_mass.npy"), mass)
    np.save(os.path.join(INPD, "_L2_p21pos.npy"), a)

    # 副法: 单步 24h→7d
    idx_7d_p16 = ((adata.obs["day"] == 7.0) & (adata.obs["Cdkn2a"] > 0)).values
    P1s = P1.copy(); P1s = sparse.csr_matrix(P1s)
    rs1 = np.asarray(P1s.sum(axis=1)).ravel(); rs1[rs1 == 0] = 1
    P1s = sparse.diags(1/rs1) @ P1s
    m7 = np.asarray(P1s[np.ix_(idx_ac, np.where(idx_7d_p16)[0])].sum(axis=1)).ravel()
    R7 = m7[a].mean() / max(m7[~a].mean(), 1e-12)
    results["moscot_step7d"] = dict(R=R7)
    log(f"  24h→7d p16+ 单步 R = {R7:.3f} (7d p16+ n={idx_7d_p16.sum()})")

except Exception as e:
    import traceback; traceback.print_exc()
    log(f"  !! moscot 失败: {e}")

# ===== 判定 =====
log("\n" + "=" * 72)
log("判定 (预登记):")
for k, r in results.items():
    R, p = r.get("R", np.nan), r.get("p_perm", np.nan)
    if np.isnan(R): continue
    v = "延续" if (R >= 2 and (not np.isnan(p)) and p < 0.05) else ("独立起源" if 0.5 <= R <= 2 else "中间")
    log(f"  [{k}] R={R:.3f}, p={p if not np.isnan(p) else 'N/A'} → {v}")
log(f"\n[L2c] DONE ({time.time()-t0:.0f}s)")
with open(os.path.join(INPD, "L2_results.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
pd.DataFrame([{"method": k, **v} for k, v in results.items()]).to_csv(
    os.path.join(INPD, "L2_enrichment.csv"), index=False)
