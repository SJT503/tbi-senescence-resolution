# s1_16b_decont_stats.py — 去污染后 GSE277487 时间轴统计 (补审计 B1 的统计端)
# 输入: results/s1_277487_decont/panel_decont_allcells.csv.gz (s1_16a 产出, decontXcounts)
# 口径: 与 s1_11b 逐字一致 —— 每型自己的共同 floor, 解析下采样 pdet, 鼠级(>=10细胞), KW + 成对 exact MWU
#       新增: p16 全型同场 + 身份基因自检 + 统一 floor 敏感性
# 用法: python s1_16b_decont_stats.py
# 产物: results/s1_277487_decont/{stats_permouse.csv, stats_summary.csv, LOG_STATS.txt}
import os
import numpy as np, pandas as pd
from scipy.special import gammaln
from scipy.stats import kruskal, mannwhitneyu

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
INP  = os.path.join(SEN, "results/s1_277487_decont/panel_decont_allcells.csv.gz")
OUTD = os.path.join(SEN, "results/s1_277487_decont")
TPS  = ["Ctrl", "group1", "group2", "group3"]
LBL  = {"Ctrl": "Ctrl", "group1": "6h", "group2": "2d", "group3": "4d"}
TYPES = ["Oligo","Astrocyte","Microglia","ExN","InN","OPC","Pericyte",
         "Ambiguous_Endothel_ExN","Ambiguous_ChoroidPlx_Pericyte"]

L = []
def log(s=""):
    print(s, flush=True); L.append(s)

def lc(n, k):
    return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)

def pdet(k, N, F):
    """P(下采样到F后至少保留1个); N<=F 保持原样 —— 与 s1_11b 同一函数"""
    k = np.asarray(k, float); N = np.asarray(N, float)
    res = np.where(N <= F, (k > 0).astype(float), np.nan)
    dn = N > F
    if dn.any():
        kk, NN = k[dn], N[dn]
        p = np.where(kk <= 0, 0.0, np.nan)
        nz = kk > 0
        if nz.any():
            p[nz] = 1.0 - np.exp(lc(NN[nz]-kk[nz], F) - lc(NN[nz], F))
        res[dn] = p
    return np.clip(np.nan_to_num(res, nan=0.0), 0, 1)

df = pd.read_csv(INP)
log(f"[s1_16b] {len(df)} 细胞 | 类型 {df.cell_type.nunique()} | 时点 {sorted(df.time_point.unique())}")

# ---- 身份自检 (吸取 A1 教训: 先证数据对再算) ----
log("\n=== 身份基因自检 (去污染后, 原始检出率) ===")
IDENT = [("Oligo","Plp1"),("Microglia","Hexb"),("Astrocyte","Aqp4"),("ExN","Snap25")]
for ty, g in IDENT:
    s = df[df.cell_type == ty]
    if len(s):
        log(f"  {ty:>10s} {g}+ = {100*(s[g]>0).mean():.1f}%  (n={len(s)})")
ok_oligo = 100*(df[df.cell_type=='Oligo'].Plp1 > 0).mean()
assert ok_oligo > 60, "身份自检失败: Oligo Plp1 <60% —— 数据有问题, 停止"

def run_analysis(floor_mode, tag):
    log(f"\n{'='*84}\n########## floor 模式: {tag} ##########")
    rows_sum, rows_pm = [], []
    for ty in TYPES:
        d0 = df[df.cell_type == ty]
        if len(d0) < 300: continue
        med = d0.groupby("time_point")["nUMI"].median()
        if not all(t in med.index for t in TPS): continue
        if floor_mode == "per_type":
            F = int(np.floor(med[TPS].min()))
        else:
            F = floor_mode  # unified
        d = d0.copy()
        d["p21"] = 100*pdet(d.Cdkn1a.values, d.nUMI.values, F)
        d["p16"] = 100*pdet(d.Cdkn2a.values, d.nUMI.values, F)
        g = d.groupby(["sample","time_point"]).agg(
                p21=("p21","mean"), p16=("p16","mean"), n=("p21","size")).reset_index()
        g = g[g.n >= 10]
        cnt_tp = g.groupby("time_point").size()
        if any(cnt_tp.get(t, 0) < 3 for t in TPS): continue
        g["celltype"] = ty
        rows_pm.append(g)
        for met in ["p21","p16"]:
            by = {t: g.loc[g.time_point==t, met].values for t in TPS}
            if np.allclose(np.concatenate([by[t] for t in TPS]), 0):
                continue
            kw = kruskal(*[by[t] for t in TPS]).pvalue
            ps = {}
            for t in TPS[1:]:
                try:
                    ps[t] = mannwhitneyu(by[t], by["Ctrl"], alternative="two-sided", method="exact").pvalue
                except Exception:
                    ps[t] = np.nan
            base = max(by["Ctrl"].mean(), 1e-9)
            rows_sum.append(dict(celltype=ty, metric=met, floor=F, F_or_base=f"Ctrl={by['Ctrl'].mean():.3f}",
                h6=by["group1"].mean(), d2=by["group2"].mean(), d4=by["group3"].mean(),
                fold6=by["group1"].mean()/base, fold4=by["group3"].mean()/base,
                diff4=by["group3"].mean()-by["Ctrl"].mean(),
                KW_p=kw, p6=ps["group1"], p4=ps["group3"],
                n_mice="/".join(str(cnt_tp[t]) for t in TPS)))
    s = pd.DataFrame(rows_sum)
    pd.concat(rows_pm, ignore_index=True).to_csv(os.path.join(OUTD, f"stats_permouse_{tag}.csv"), index=False)
    if len(s):
        s.to_csv(os.path.join(OUTD, f"stats_summary_{tag}.csv"), index=False)
        for met in ["p21","p16"]:
            sub = s[s.metric==met].sort_values("fold6" if met=="p21" else "diff4", ascending=False)
            if not len(sub): continue
            log(f"\n--- {met} (按{'6h倍数' if met=='p21' else '4d绝对差'}降序) ---")
            log(f"{'celltype':>28s} {'floor':>6s} {'Ctrl':>7s} {'6h':>8s} {'2d':>8s} {'4d':>8s} {'6h倍':>6s} {'4d倍':>6s} {'4d-Ctrl':>8s} {'KW_p':>8s} {'p6':>7s} {'p4':>7s}")
            for _, r in sub.iterrows():
                ctrl = float(r.F_or_base.split("=")[1])
                log(f"{r.celltype:>28s} {r.floor:6d} {ctrl:7.3f} {r.h6:8.3f} {r.d2:8.3f} {r.d4:8.3f} "
                    f"{r.fold6:6.2f} {r.fold4:6.2f} {r.diff4:+8.3f} {r.KW_p:8.3g} {r.p6:7.3g} {r.p4:7.3g}")
    return s

# 每型自己的 floor (与 s1_11b 可比)
s1 = run_analysis("per_type", "pertype")
# 统一 floor (最小者) —— 小胶质持续性的稳健性检验
meds = []
for ty in TYPES:
    d0 = df[df.cell_type==ty]
    if len(d0) < 300: continue
    m = d0.groupby("time_point")["nUMI"].median()
    if all(t in m.index for t in TPS): meds.append(int(np.floor(m[TPS].min())))
UF = min(meds)
log(f"\n统一 floor = {UF} (各型: {sorted(meds)})")
s2 = run_analysis(UF, "unified")

log("\n" + "="*84)
log("关键判读 (对照旧 raw 版 s1_11b: Oligo 6h 20.46/7.00x; Microglia 6h 8.62/36.6x 4d 5.55x; Astro 12.4x)")
for tag, s in [("pertype", s1), ("unified", s2)]:
    for ty in ["Oligo","Microglia","Astrocyte"]:
        r = s[(s.celltype==ty)&(s.metric=="p21")]
        if len(r):
            r = r.iloc[0]
            log(f"  [{tag}] {ty:>10s}: 6h {r.h6:6.2f} ({r.fold6:5.2f}x, p={r.p6:.3g}) | 4d {r.d4:5.2f} ({r.fold4:5.2f}x) | 4d-Ctrl {r.diff4:+.3f}")
mg16 = s1[(s1.celltype=="Microglia")&(s1.metric=="p16")]
if len(mg16):
    r = mg16.iloc[0]
    log(f"  [p16/小胶质 去污染后] Ctrl {float(r.F_or_base.split('=')[1]):.3f} | 6h {r.h6:.3f} | 4d {r.d4:.3f} (p6={r.p6:.3g}, p4={r.p4:.3g})")

with open(os.path.join(OUTD, "LOG_STATS.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
log(f"\n[s1_16b] DONE → {OUTD}")
