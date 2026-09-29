# s1_11b_celltype_stats.py — 跨细胞型 p21 时间轴统计 (纯解析下采样, 不抽样)
# 输入: results/s1_277487_timecourse/p21_panel_allcells.csv.gz (s1_11a 导出, 读一次 ckpt)
# 方法: 解析下采样概率 —— 细胞原有 k 个 p21 分子、N 个 UMI, 下采样到 F 个 UMI 时
#        P(至少抽到 1 个) = 1 - C(N-k, F)/C(N, F)  (超几何)
#        这是"无限次抽样的期望检出率", 比单次随机实现更准且无随机性
# 统计: 每型共同 floor = 该型四时点中位 nUMI 的最小值; 鼠级(每鼠>=10细胞) 检出率
#       非单调且四水平 → Kruskal-Wallis; 各时点 vs Ctrl → Wilcoxon exact; BH 校正
# 用法: python s1_11b_celltype_stats.py
# 产物: results/s1_277487_timecourse/{celltype_p21_summary.csv, celltype_p21_permouse.csv, LOG_CELLTYPE.txt}
import os, math, sys
import numpy as np, pandas as pd
from scipy.stats import kruskal, mannwhitneyu, false_discovery_control

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
INP  = os.path.join(SEN, "results/s1_277487_timecourse/p21_panel_allcells.csv.gz")
OUTD = os.path.join(SEN, "results/s1_277487_timecourse")

TPS  = ["Ctrl", "group1", "group2", "group3"]
LBL  = {"Ctrl": "Ctrl", "group1": "6h", "group2": "2d", "group3": "4d"}
MIN_CELLS_PER_MOUSE = 10
MIN_MICE = 3
TYPES = ["Oligo","Astrocyte","Microglia","ExN","InN","OPC","Pericyte",
         "Ependymal","Proliferat","Ambiguous_Endothel_ExN","Ambiguous_ChoroidPlx_Pericyte"]

L = []
def log(s=""):
    print(s, flush=True); L.append(s)

# ---------- 解析下采样检出概率 (log 空间, 防组合数上溢) ----------
def log_comb(n, k):
    """log C(n,k), 支持数组 k; n 标量或数组"""
    n = np.asarray(n, dtype=float); k = np.asarray(k, dtype=float)
    # C(n,k) = exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1))
    from scipy.special import gammaln
    out = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    return out

def p_detect_after_downsample(k, N, F):
    """P(至少抽到1个) = 1 - C(N-k,F)/C(N,F); 若 N<=F 则必保留 → 1 if k>0 else 0"""
    k = np.asarray(k, dtype=float); N = np.asarray(N, dtype=float)
    F = np.asarray(F, dtype=float)
    res = np.zeros_like(N, dtype=float)
    # 未下采样 (N<=F): 原始检出
    nodown = N <= F
    res[nodown] = (k[nodown] > 0).astype(float)
    # 下采样: N>F
    dn = ~nodown
    if np.any(dn):
        kk, NN = k[dn], N[dn]
        # k==0 → 必检不到
        zero = kk <= 0
        p = np.zeros_like(NN)
        nz = ~zero
        if np.any(nz):
            # 正确形式: p = 1 - C(N-k,F)/C(N,F)  →  对数比 = logC(N-k,F) - logC(N,F)
            lg = log_comb(NN[nz] - kk[nz], F) - log_comb(NN[nz], F)
            p[nz] = 1.0 - np.exp(lg)      # C(N-k,F)/C(N,F)
        res[dn] = p
    return np.clip(res, 0.0, 1.0)

# ---------- 读入小表 ----------
log(f"[s1_11b] 读入 {INP}")
df = pd.read_csv(INP)
log(f"[s1_11b] 全部细胞 {len(df)} | 类型 {df.cell_type.nunique()} | 时点 {sorted(df.time_point.unique())}")
log(f"           全库 Cdkn1a 阳性率 = {100*(df.Cdkn1a>0).mean():.2f}%")

rows_sum, rows_pm = [], []
for ty in TYPES:
    d = df[df.cell_type == ty].copy()
    if len(d) < 50:
        log(f"  skip {ty:32s} cells={len(d)}")
        continue
    med = d.groupby("time_point")["nUMI"].median()
    if not all(t in med.index for t in TPS):
        log(f"  skip {ty:32s} (缺时点)")
        continue
    F = int(np.floor(med[TPS].min()))
    d["p21_det"] = 100.0 * p_detect_after_downsample(d.Cdkn1a.values, d.nUMI.values, F)
    d["p16_det"] = 100.0 * p_detect_after_downsample(d.Cdkn2a.values, d.nUMI.values, F)
    # 鼠级 (>=MIN_CELLS_PER_MOUSE)
    g = d.groupby(["sample", "time_point"]).agg(
            p21=("p21_det", "mean"), p16=("p16_det", "mean"), n=("p21_det", "size")).reset_index()
    g = g[g.n >= MIN_CELLS_PER_MOUSE]
    cnt_tp = g.groupby("time_point").size()
    if any(cnt_tp.get(t, 0) < MIN_MICE for t in TPS):
        log(f"  skip {ty:32s} (合格鼠不足: {dict(cnt_tp)})")
        continue
    g["celltype"] = ty
    rows_pm.append(g)
    by = {t: g.loc[g.time_point == t, "p21"].values for t in TPS}
    kw = kruskal(*[by[t] for t in TPS]).pvalue
    ps = {}
    for t in TPS[1:]:
        try:
            ps[t] = mannwhitneyu(by[t], by["Ctrl"], alternative="two-sided", method="exact").pvalue
        except Exception:
            ps[t] = np.nan
    rows_sum.append(dict(celltype=ty, floor=F, n_mice="/".join(str(cnt_tp[t]) for t in TPS),
                         Ctrl=np.mean(by["Ctrl"]), h6=np.mean(by["group1"]),
                         d2=np.mean(by["group2"]), d4=np.mean(by["group3"]),
                         fold_6h=np.mean(by["group1"]) / max(np.mean(by["Ctrl"]), 1e-9),
                         KW_p=kw, p_6h=ps["group1"], p_2d=ps["group2"], p_4d=ps["group3"]))
    log(f"  {ty:32s} F={F:5d} | Ctrl {np.mean(by['Ctrl']):7.3f} | 6h {np.mean(by['group1']):7.3f} "
        f"({rows_sum[-1]['fold_6h']:5.2f}x) | 2d {np.mean(by['group2']):7.3f} | 4d {np.mean(by['group4' if False else 'group3']):7.3f} | KW p={kw:.4g}")

sm = pd.DataFrame(rows_sum)
if len(sm) == 0:
    log("[s1_11b] 没有合格的细胞型, 退出"); sys.exit(1)
# BH 校正 6h 检验
sm["p_6h_BH"] = false_discovery_control(sm.p_6h.fillna(1.0).values)
sm = sm.sort_values("fold_6h", ascending=False).reset_index(drop=True)
sm.to_csv(os.path.join(OUTD, "celltype_p21_summary.csv"), index=False)
pd.concat(rows_pm, ignore_index=True).to_csv(os.path.join(OUTD, "celltype_p21_permouse.csv"), index=False)

log("\n" + "=" * 84)
log("汇总 (按 6h 倍数降序):")
log(sm[["celltype","Ctrl","h6","d2","d4","fold_6h","KW_p","p_6h","p_6h_BH"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))

ol = sm[sm.celltype == "Oligo"]
if len(ol):
    rank = int((sm.fold_6h > ol.fold_6h.values[0]).sum()) + 1
    oth = sm[sm.celltype != "Oligo"].fold_6h
    log(f"\n关键判读:")
    log(f"  Oligo fold_6h = {ol.fold_6h.values[0]:.2f}  →  排名 {rank}/{len(sm)}")
    log(f"  其他型 fold_6h 中位 = {oth.median():.2f}  (范围 {oth.min():.2f}–{oth.max():.2f})")
    log(f"  其他型中 fold_6h >= 2 的个数 = {(oth >= 2).sum()} / {len(oth)}")
    log(f"  → 少突特异性判定: {'支持(Oligo 第一且其余接近1)' if rank == 1 else '需谨慎(非唯一第一)'}")
else:
    log("\n!! Oligo 未进入汇总表, 检查入组条件")

with open(os.path.join(OUTD, "LOG_CELLTYPE.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
log(f"\n[s1_11b] DONE → {OUTD}")
