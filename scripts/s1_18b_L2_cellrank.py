# -*- coding: utf-8 -*-
"""s1_18b_L2_cellrank.py — L2 起源裁决本体 (预登记: L2_PROTOCOL.md)
主法: CellRank2 RealTimeKernel (真实时点 OT) — 富集比 R = P(6mo p16+ | 急性 p21+) / P(6mo p16+ | 急性 p21−)
副法: moscot TemporalProblem 链式传输 (24h→7d→6mo)
判据: R>=2 & p_perm<0.05 → 延续; R in 0.5-2 & ns → 独立起源
产物: results/s1_L2_origin/{L2_results.txt, L2_enrichment.csv}
"""
import os, sys, time
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
CAP_PER_TP = 6000   # 计算必需的截断(OT CPU), 固定seed, 记录在输出; 判据不变

t0 = time.time()
log("[L2] 读入矩阵 ...")
X = mmread(os.path.join(INPD, "counts.mtx")).T.tocsr()          # cells x genes
genes = [l.strip() for l in open(os.path.join(INPD, "genes.txt"))]
cells = pd.read_csv(os.path.join(INPD, "cells.csv"), index_col="cell_id")
assert X.shape == (len(cells), len(genes)), f"{X.shape} vs {(len(cells), len(genes))}"
adata = ad.AnnData(X=X, obs=cells, var=pd.DataFrame(index=genes))
adata.obs["day"] = pd.Categorical(adata.obs["time_point"].map(DAY), ordered=True)  # RealTimeKernel 要求 categorical
log(f"[L2] {adata.n_obs} 细胞 x {adata.n_vars} 基因 | {dict(adata.obs.time_point.value_counts())}")

# 计算截断 (每时点 <=6000)
keep = []
for tp, idx in adata.obs.groupby("time_point").groups.items():
    idx = np.asarray(idx)
    if len(idx) > CAP_PER_TP:
        idx = np.random.choice(idx, CAP_PER_TP, replace=False)
        log(f"  截断 {tp}: {len(adata.obs.loc[adata.obs.time_point==tp])} → {CAP_PER_TP}")
    keep += list(idx)
adata = adata[keep].copy()

# ---- 预处理 ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, svd_solver="arpack", random_state=SEED)
try:
    sc.external.pp.harmony_integrate(adata, key="lib", basis="X_pca",
                                     adjusted_basis="X_pca_harmony", max_iter_harmony=20)
    rep = "X_pca_harmony"; log("[L2] Harmony 校正完成 (lib)")
except Exception as e:
    rep = "X_pca"; log(f"[L2] Harmony 失败({e}), 用原始 PCA — 跨文库残差风险声明")
sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, random_state=SEED)
log(f"[L2] 预处理完成 ({time.time()-t0:.0f}s)")

# 状态旗标 (清洗计数检出式, 来自 s1_18a cells.csv; 截断后仍在 obs)
acute = (adata.obs["day"] == 1.0).values
term_p16 = ((adata.obs["day"] == 180.0) & (adata.obs["Cdkn2a"] > 0)).values
p21pos = (adata.obs["Cdkn1a"] > 0).values
log(f"[L2] 急性(24h) n={acute.sum()} | p21+ {p21pos[acute].sum()} | 6mo p16+ n={term_p16.sum()}")

def enrichment_from_T(T, acute, p21pos, term_p16):
    """R = mean(行和到 6mo p16+ | 急性 p21+) / mean(... | 急性 p21−)"""
    mass = np.asarray(T[acute][:, term_p16].sum(axis=1)).ravel()
    a = p21pos[acute]
    m_pos, m_neg = mass[a].mean(), mass[~a].mean()
    return (m_pos / m_neg if m_neg > 0 else np.nan), m_pos, m_neg

def permutation_p(T, acute, p21pos, term_p16, n=1000):
    mass = np.asarray(T[acute][:, term_p16].sum(axis=1)).ravel()
    a = p21pos[acute]
    obs_R = mass[a].mean() / max(mass[~a].mean(), 1e-12)
    rng = np.random.default_rng(SEED)
    cnt = 0
    for _ in range(n):
        perm = rng.permutation(a)
        r = mass[perm].mean() / max(mass[~perm].mean(), 1e-12)
        if r >= obs_R: cnt += 1
    return obs_R, (cnt + 1) / (n + 1)

results = {}

# ================= 主法: CellRank RealTimeKernel =================
log("\n[L2] === 主法: CellRank RealTimeKernel ===")
try:
    from cellrank.kernels import RealTimeKernel
    rt = RealTimeKernel(adata, time_key="day")
    rt.compute_transition_matrix()
    T = rt.transition_matrix.tocsr()
    R, mpos, mneg = enrichment_from_T(T, acute, p21pos, term_p16)
    R_perm, p_perm = permutation_p(T, acute, p21pos, term_p16)
    results["RTK"] = dict(R=R, R_perm=R_perm, p_perm=p_perm, m_pos=mpos, m_neg=mneg)
    log(f"  急性p21+→6mo p16+ 质量: {mpos:.4f} vs p21−: {mneg:.4f}")
    log(f"  富集比 R = {R:.3f} | 置换 p = {p_perm:.4f} (n=1000)")
except Exception as e:
    import traceback; traceback.print_exc()
    log(f"  !! RealTimeKernel 失败: {e}")

# ================= 副法: moscot 链式 =================
log("\n[L2] === 副法: moscot TemporalProblem (24h→7d→6mo 链式) ===")
try:
    import moscot as mc
    try:
        from moscot.problems.time_series import TemporalProblem   # moscot <=0.4
    except ImportError:
        from moscot.problems import TemporalProblem               # moscot >=0.5
    tp = TemporalProblem(adata)
    tp = tp.prepare(time_key="day")
    tp = tp.solve(epsilon=0.005, tau_a=0.99, tau_b=0.99, seed=SEED)
    P1 = tp[1.0, 7.0].solution.transport_matrix     # 24h x 7d
    P2 = tp[7.0, 180.0].solution.transport_matrix    # 7d x 6mo
    Pch = (P1 @ P2).tocsr()                          # 24h x 6mo 链式
    # 行归一 (质量守恒近似)
    rs = np.asarray(Pch.sum(axis=1)).ravel(); rs[rs == 0] = 1
    Pch = sparse.diags(1/rs) @ Pch
    idx_acute = np.where(acute)[0]
    idx_term = np.where(term_p16)[0]
    mass = np.asarray(Pch[np.ix_(idx_acute, idx_term)].sum(axis=1)).ravel()
    a = p21pos[acute]
    Rm = mass[a].mean() / max(mass[~a].mean(), 1e-12)
    rng = np.random.default_rng(SEED); cnt = 0
    for _ in range(1000):
        perm = rng.permutation(a)
        r = mass[perm].mean() / max(mass[~perm].mean(), 1e-12)
        if r >= Rm: cnt += 1
    results["moscot"] = dict(R=Rm, p_perm=(cnt+1)/1001)
    log(f"  链式富集比 R = {Rm:.3f} | 置换 p = {(cnt+1)/1001:.4f}")
except Exception as e:
    import traceback; traceback.print_exc()
    log(f"  !! moscot 失败 (回退按协议: 仅主法) : {e}")

# ================= 深度敏感性 (窗口内重算, 仅主法) =================
log("\n[L2] === 深度敏感性 (nUMI IQR 窗口, 主法重算) ===")
try:
    nUMI = np.asarray(adata.X.sum(axis=1)).ravel() if adata.X is not None else None
    # X 已 scale — 用 obs 里没有 nUMI; 用 Mcl1+Cdkn1a 等计数和作代理不可靠 → 用 neighbors 图上重算不可行
    # 简化: 用原始计数文件重读一次拿 nUMI
    Xraw = mmread(os.path.join(INPD, "counts.mtx")).T.tocsr()
    allcells = pd.read_csv(os.path.join(INPD, "cells.csv"))
    nUMI_all = np.asarray(Xraw.sum(axis=1)).ravel()
    pos_map = {c: i for i, c in enumerate(allcells.cell_id)}
    idx_orig = np.array([pos_map[c] for c in adata.obs_names])
    nUMI = nUMI_all[idx_orig]
    q1, q3 = np.percentile(nUMI, [25, 75])
    win = (nUMI >= q1) & (nUMI <= q3) & acute
    if "RTK" not in results: raise RuntimeError("主法未成功, 敏感性依赖 rt 跳过")
    Tt = rt.transition_matrix.tocsr()
    mass = np.asarray(Tt[win][:, term_p16].sum(axis=1)).ravel()
    a = p21pos[win]
    if a.sum() >= 20 and (~a).sum() >= 20:
        Rw = mass[a].mean() / max(mass[~a].mean(), 1e-12)
        log(f"  窗口内 R = {Rw:.3f} (n+={a.sum()}, n−={(~a).sum()})")
        results.setdefault("RTK", {})["R_window"] = Rw
except Exception as e:
    log(f"  敏感性失败: {e}")

# ================= 判定 (按 L2_PROTOCOL 预登记) =================
log("\n" + "=" * 72)
log("判定 (预登记判据):")
for k, r in results.items():
    if "R" not in r: continue
    R, p = r.get("R", np.nan), r.get("p_perm", np.nan)
    verdict = "延续" if (R >= 2 and p < 0.05) else ("独立起源" if (0.5 <= R <= 2) else "中间/不可判")
    log(f"  [{k}] R={R:.3f}, p_perm={p:.4f} → {verdict}")
log(f"\n[L2] DONE ({time.time()-t0:.0f}s) → {INPD}")
with open(os.path.join(INPD, "L2_results.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
pd.DataFrame([{"method": k, **v} for k, v in results.items()]).to_csv(
    os.path.join(INPD, "L2_enrichment.csv"), index=False)
