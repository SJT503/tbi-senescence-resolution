# s1_12_crosscelltype_24h.py — 24h 批跨细胞型对照: 剂量梯度 + 同侧主导
# 目的: "剂量梯度"与"同侧主导"这两个 s1_05b 只在 Oligo 上做过的检验, 在其他细胞型是否也成立
#       → 决定"少突特异性"主张能否成立(跨型分析已否掉"6h 峰特异", 这是最后一炮)
# 数据: processed_data/s1_24h_pilot/_covary_export.csv.gz (240,498 细胞, 原始 nUMI + Cdkn1a + cell_type)
# 方法: 解析下采样 (超几何闭式解, 与 s1_11b 同一函数, 已用 scipy 对拍 8/8 通过)
#       与 s1_03 的 prop 法在 N<=F 时完全一致(都=原样保留); N>F 时本式用超几何而非二项薄化
# 用法: python s1_12_crosscelltype_24h.py
# 产物: results/s1_crosscelltype_24h/{dose_gradient.csv, ipsilateral.csv, LOG.txt}
import os, math
import numpy as np, pandas as pd
from scipy.special import gammaln
from scipy.stats import mannwhitneyu, wilcoxon, fisher_exact
from statistics import NormalDist

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
INP  = os.path.join(SEN, "processed_data/s1_24h_pilot/_covary_export.csv.gz")
OUTD = os.path.join(SEN, "results/s1_crosscelltype_24h"); os.makedirs(OUTD, exist_ok=True)
FLOOR = 1921          # 与 s1_03 协议一致
MIN_CELLS = 10

L = []
def log(s=""):
    print(s, flush=True); L.append(s)

def log_comb(n, k):
    return gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)

def p_detect(k, N, F):
    """P(下采样到 F 后至少保留 1 个) = 1 - C(N-k,F)/C(N,F); N<=F 时原样(检出 iff k>0)
    —— 与 s1_03 的 prop=ifelse(tot>F, F/tot, 1) 在 N<=F 分支完全一致"""
    k = np.asarray(k, float); N = np.asarray(N, float); F = float(F)
    res = np.zeros_like(N)
    nod = N <= F
    res[nod] = (k[nod] > 0).astype(float)
    dn = ~nod
    if np.any(dn):
        kk, NN = k[dn], N[dn]
        p = np.zeros_like(NN)
        nz = kk > 0
        if np.any(nz):
            lg = log_comb(NN[nz] - kk[nz], F) - log_comb(NN[nz], F)
            p[nz] = 1.0 - np.exp(lg)
        res[dn] = p
    return np.clip(res, 0.0, 1.0)

def cliff(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return (sum((x > b).sum() for x in a) - sum((x < b).sum() for x in a)) / (len(a) * len(b))

def jt_increasing(groups_vals):
    """Jonckheere 递增趋势: J<EJ → z<0 → 单侧 p=Φ(z) (沿用 s1_05b 修正版)"""
    ns = [len(v) for v in groups_vals]
    J = sum(((a[:, None] > b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum())
            for i, a in enumerate(groups_vals) for b in groups_vals[i+1:])
    N = sum(ns)
    EJ = (N**2 - sum(n*n for n in ns)) / 4
    VJ = (N**2*(2*N+3) - sum(n*n*(2*n+3) for n in ns)) / 72
    z = (J - EJ) / math.sqrt(VJ) if VJ > 0 else np.nan
    return J, z, NormalDist().cdf(z)

# ---- 分组 (与 s1_05b 逐字一致) ----
GROUPS = dict(
    naive_M=["RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25"],
    CCI_ipsi_M=["RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L"],
    CCI_HS_M=["RNASEQ01","RNASEQ07","RNASEQ13"],
    rCHI_M=["RNASEQ09","RNASEQ11","RNASEQ12"],
    CCI_contra_M=["RNASEQ17R","RNASEQ18R","RNASEQ19R"],
    CCI_F=["RNASEQ20L","RNASEQ21L","RNASEQ22L"],
    naive_F=["RNASEQ26","RNASEQ27","RNASEQ28"])
SID2GRP = {s: g for g, v in GROUPS.items() for s in v}
# s1_05b S1 用的子批池 (保持逐位一致)
DOSE = [["RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16"],   # naive 子批A
        ["RNASEQ09","RNASEQ11","RNASEQ12"],               # rCHI
        ["RNASEQ03","RNASEQ06","RNASEQ15"],               # CCI 子批A
        ["RNASEQ01","RNASEQ07","RNASEQ13"]]               # CCI+HS
PAIRS = [("RNASEQ17L","RNASEQ17R"), ("RNASEQ18L","RNASEQ18R"), ("RNASEQ19L","RNASEQ19R")]

df = pd.read_csv(INP)
log(f"[s1_12] 细胞 {len(df)} | 类型 {df.cell_type.nunique()} | floor={FLOOR}")
df["grp"] = df["sample"].map(SID2GRP)
df["det"] = 100.0 * p_detect(df.Cdkn1a.values, df.nUMI.values, FLOOR)

# 逐鼠逐型检出率
pm = (df.groupby(["cell_type", "sample", "grp"])
        .agg(n=("det", "size"), det=("det", "mean")).reset_index())
pm = pm[pm.n >= MIN_CELLS]
pm.to_csv(os.path.join(OUTD, "permouse_alltypes.csv"), index=False)

# ---- 检验 1: 剂量梯度 (naive < rCHI < CCI_ipsi < CCI_HS) ----
log("\n" + "=" * 78)
log("检验1 | 剂量梯度 (naive_A < rCHI < CCI_A < CCI+HS), Jonckheere 单侧递增")
log(f"{'cell_type':>26s} {'n梯度':>5s} {'naive':>8s} {'rCHI':>8s} {'CCI':>8s} {'CCI+HS':>8s} {'J':>5s} {'z':>6s} {'p':>8s}")
rows1 = []
for ty in sorted(pm.cell_type.unique()):
    seq = []
    for sids in DOSE:
        v = pm[(pm.cell_type == ty) & (pm["sample"].isin(sids))]["det"].values
        seq.append(v)
    if any(len(v) == 0 for v in seq):
        continue
    J, z, p = jt_increasing(seq)
    means = [v.mean() for v in seq]
    rows1.append(dict(cell_type=ty, n_per_group="/".join(str(len(v)) for v in seq),
                      naive=means[0], rCHI=means[1], CCI=means[2], CCI_HS=means[3],
                      J=J, z=z, p_trend=p))
    log(f"{ty:>26s} {'/'.join(str(len(v)) for v in seq):>5s} "
        f"{means[0]:8.3f} {means[1]:8.3f} {means[2]:8.3f} {means[3]:8.3f} {J:5.0f} {z:6.2f} {p:8.4g}")
d1 = pd.DataFrame(rows1)
if len(d1):
    d1 = d1.sort_values("p_trend")
    d1.to_csv(os.path.join(OUTD, "dose_gradient.csv"), index=False)

# ---- 检验 2: 同侧主导 (配对 ipsi vs contra) ----
log("\n" + "=" * 78)
log("检验2 | 同侧主导 (同鼠配对 17L/18L/19L vs 17R/18R/19R)")
log(f"{'cell_type':>26s} {'可用对':>6s} {'ipsi':>9s} {'contra':>9s} {'升幅%':>9s} {'Wilcoxon p':>11s}")
rows2 = []
for ty in sorted(pm.cell_type.unique()):
    pairs = []
    for L_, R_ in PAIRS:
        a = pm[(pm.cell_type == ty) & (pm["sample"] == L_)]["det"].values
        b = pm[(pm.cell_type == ty) & (pm["sample"] == R_)]["det"].values
        if len(a) and len(b):
            pairs.append((a[0], b[0]))
    if len(pairs) < 2:
        rows2.append(dict(cell_type=ty, n_pairs=len(pairs), note="不可评(<2对)"))
        log(f"{ty:>26s} {len(pairs):>6d} {'':>9s} {'':>9s} {'':>9s} {'不可评(<2对)':>11s}")
        continue
    ip = [x for x, _ in pairs]; ct = [y for _, y in pairs]
    try:
        pw = wilcoxon(ip, ct).pvalue if len(ip) >= 3 else float("nan")
    except ValueError:
        pw = float("nan")
    lift = [100 * (x - y) / y if y > 0 else np.nan for x, y in pairs]
    rows2.append(dict(cell_type=ty, n_pairs=len(pairs), ipsi=np.mean(ip), contra=np.mean(ct),
                      lift_pct=np.nanmean(lift), wilcoxon_p=pw))
    log(f"{ty:>26s} {len(pairs):>6d} {np.mean(ip):9.3f} {np.mean(ct):9.3f} "
        f"{np.nanmean(lift):9.1f} {pw:11.4g}")
d2 = pd.DataFrame(rows2)
d2.to_csv(os.path.join(OUTD, "ipsilateral.csv"), index=False)

# ---- 关键判读 ----
log("\n" + "=" * 78)
log("关键判读:")
if len(d1):
    ol = d1[d1.cell_type == "Oligo"]
    if len(ol):
        r = int((d1.p_trend < ol.p_trend.values[0]).sum()) + 1
        log(f"  剂量梯度: Oligo p={ol.p_trend.values[0]:.4g}, 排名 {r}/{len(d1)} (越小越显著)")
        sig = d1[d1.p_trend < 0.05]
        log(f"  剂量梯度显著(p<0.05)的型: {', '.join(sig.cell_type) if len(sig) else '无'}")
if len(d2) and "wilcoxon_p" in d2.columns:
    ol2 = d2[d2.cell_type == "Oligo"]
    if len(ol2) and pd.notna(ol2.lift_pct.values[0]):
        r2 = int((d2.lift_pct > ol2.lift_pct.values[0]).sum()) + 1
        log(f"  同侧升幅: Oligo {ol2.lift_pct.values[0]:.1f}%, 排名 {r2}/{d2.n_pairs.notna().sum()}")
log(f"\n  [对照] 上一轮 s1_05b 在 Oligo 上: 剂量梯度 p=0.00503 | 同侧升幅 +424%")

with open(os.path.join(OUTD, "LOG.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
log(f"\n[s1_12] DONE → {OUTD}")
