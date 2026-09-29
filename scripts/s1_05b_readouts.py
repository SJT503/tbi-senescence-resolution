# -*- coding: utf-8 -*-
"""s1_05b_readouts.py — Step-1 第5步b: 预登记读出 P1-P3 / S1-S3 (协议§4)

口径逐位对齐 planA_47_audit_gapfill.py (G7 bg7 / G8 Fisher OR / G4 exact MWU+Cliff)。
输入: processed_data/s1_24h_pilot/{readout_celltable.csv.gz, labels/all_cells_labels.csv.gz}
输出: results/s1_24h_pilot/readout_summary.txt + readout_permouse.csv + readout_p2bytype.csv
用法: python s1_05b_readouts.py [version]   (默认 main; lo20/hi20 为敏感性)
"""
import os, sys, math, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon, fisher_exact

SEN = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
CT  = os.path.join(SEN, "processed_data/s1_24h_pilot/readout_celltable.csv.gz")
LB  = os.path.join(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")
OUTD = os.path.join(SEN, "results/s1_24h_pilot"); os.makedirs(OUTD, exist_ok=True)

VER = sys.argv[1] if len(sys.argv) > 1 else "main"
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

def cliff(a, b):  # 逐字同 planA_47
    a, b = np.asarray(a, float), np.asarray(b, float)
    return (sum((x > b).sum() for x in a) - sum((x < b).sum() for x in a)) / (len(a) * len(b))

def mwu(x, y, tag):
    p = mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue
    log(f"    {tag}: inj mean={np.mean(x):.4g} (n={len(x)}) vs sham mean={np.mean(y):.4g} (n={len(y)})"
        f"  exact p={p:.4g}  δ={cliff(np.asarray(x), np.asarray(y)):+.2f}")
    return p

GROUPS = dict(
    naive_M=["RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25"],
    CCI_ipsi_M=["RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L"],
    CCI_HS_M=["RNASEQ01","RNASEQ07","RNASEQ13"],
    rCHI_M=["RNASEQ09","RNASEQ11","RNASEQ12"],
    CCI_contra_M=["RNASEQ17R","RNASEQ18R","RNASEQ19R"],
    CCI_F=["RNASEQ20L","RNASEQ21L","RNASEQ22L"],
    naive_F=["RNASEQ26","RNASEQ27","RNASEQ28"])
DIRECT13 = GROUPS["naive_M"] + GROUPS["CCI_ipsi_M"]
BG7 = ["Astrocyte","InN","Endothel","Pericyte","ChoroidPlx","Ependymal","ExN"]  # 焊死: 含ExN不含Tcell

ct = pd.read_csv(CT); ct = ct[ct.version == VER]
lb = pd.read_csv(LB, usecols=["cell_id","cell_type"])
df = ct.merge(lb, on="cell_id", how="left")
df = df[df.cell_type.notna() & (df.cell_type != "Ambiguous")]
log(f"[s1_05b] version={VER}: {len(df)} cells, {df.cell_type.nunique()} types (Ambiguous 已剔)")

MIN_CELLS = 10
df["det_p21"] = df.p21_umis.fillna(0) > 0
df["det_bax"] = df.bax_umis.fillna(0) > 0
df["det_mcl1"] = df.mcl1_umis.fillna(0) > 0
pm = (df.groupby(["sample","cell_type"])
        .agg(n=("senmayo","size"), senmayo=("senmayo","mean"),
             pct_p21=("det_p21","mean"), pct_bax=("det_bax","mean"), pct_mcl1=("det_mcl1","mean"))
        .reset_index())
pm = pm[pm.n >= MIN_CELLS]
pm["group"] = pm["sample"].map({s:g for g,v in GROUPS.items() for s in v})
pm.to_csv(os.path.join(OUTD, f"readout_permouse_{VER}.csv"), index=False)

g = lambda grp, typ, col: pm[(pm.group==grp)&(pm.cell_type==typ)][col].values
has = lambda grp, typ: len(pm[(pm.group==grp)&(pm.cell_type==typ)]) > 0

# ================= P1: oligo p21 检出率 + SenMayo, CCI vs naive =================
log("="*76); log(f"P1 | oligo Cdkn1a 检出率与 SenMayo (CCI_ipsi_M vs naive_M) [ver={VER}]")
assert has("CCI_ipsi_M","Oligo") and has("naive_M","Oligo"), "P1: Oligo 每鼠≥10细胞不满足"
p1 = {}
p1["p21"] = mwu(g("CCI_ipsi_M","Oligo","pct_p21"), g("naive_M","Oligo","pct_p21"), "pct_p21+")
p1["sen"] = mwu(g("CCI_ipsi_M","Oligo","senmayo"), g("naive_M","Oligo","senmayo"), "SenMayo_real")

log("  检出偏倚标尺 (G8 Fisher OR 复算链: 池化细胞, inj vs sham):")
ol = df[df.cell_type=="Oligo"]
inj = ol[ol.group=="CCI_ipsi_M"]; shm = ol[ol.group=="naive_M"]
ors = {}
for gene, col in [("Cdkn1a","det_p21"), ("Bax","det_bax"), ("Mcl1","det_mcl1")]:
    a, b = int(inj[col].sum()), len(inj)-int(inj[col].sum())
    c, d = int(shm[col].sum()), len(shm)-int(shm[col].sum())
    orv, pf = fisher_exact([[a,b],[c,d]])
    ors[gene] = orv
    log(f"    {gene:7s}: {a}/{len(inj)} ({100*a/len(inj):.2f}%) vs {c}/{len(shm)} ({100*c/len(shm):.2f}%)"
        f"  OR={orv:.2f}  Fisher p={pf:.3g}")
if "Cdkn1a" in ors and ors.get("Bax") is not None and ors["Bax"] > 0:
    log(f"    p21 特异性: OR_p21/OR_Bax = {ors['Cdkn1a']/max(ors['Bax'],1e-9):.1f}× "
        f"| OR_p21/OR_Mcl1 = {ors['Cdkn1a']/max(ors.get('Mcl1',float('nan')),1e-9):.1f}×  (预登记门槛 ≥3×)")

# ================= P2: 全细胞型 SenMayo 升高计数 =================
log("="*76); log(f"P2 | 各型 SenMayo 24h vs naive (鼠级 exact MWU) [ver={VER}]")
rows = []
for typ in sorted(pm.cell_type.unique()):
    x, y = g("CCI_ipsi_M",typ,"senmayo"), g("naive_M",typ,"senmayo")
    if len(x) < 3 or len(y) < 3:
        rows.append([typ,len(x),len(y),np.nan,np.nan,"n不足(<3鼠)"]); continue
    p = mannwhitneyu(x, y, alternative="greater", method="exact").pvalue
    rows.append([typ,len(x),len(y),np.mean(x)-np.mean(y),p,""])
p2 = pd.DataFrame(rows, columns=["cell_type","n_inj","n_shm","diff_mean","p_greater","note"])
p2["sig"] = p2.p_greater < 0.05
clean = p2[~p2.cell_type.str.startswith(("LowConf","Ambiguous"))]
log(f"  全部合格型 {len(p2)}: 显著升高 {int(p2.sig.sum())} | "
    f"非LowConf/Ambiguous型 {len(clean)}: 显著升高 {int(clean.sig.sum())}  (预登记问: ≥8/13?)")
log(p2.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
p2.to_csv(os.path.join(OUTD, f"readout_p2bytype_{VER}.csv"), index=False)

# ================= P3: bg7 富集 (G7 逐位复刻, 仅13直连样本) =================
log("="*76); log(f"P3 | log2(oligo/bg7) 24h vs Ctrl — 仅13直连样本 [ver={VER}]  (旧值 p=0.0734)")
p3 = pm[pm["sample"].isin(DIRECT13)]
w = p3.pivot_table(index="sample", columns="cell_type", values="senmayo")
tpv = pd.Series({s: "24h" if s in GROUPS["CCI_ipsi_M"] else "Ctrl" for s in w.index})
bg = w[[t for t in BG7 if t in w.columns]].mean(1)
rel = np.log2(w["Oligo"] / bg)
p3p = mwu(rel[tpv == "24h"].values, rel[tpv == "Ctrl"].values, "log2(oligo/bg7)")
log(f"    24h mean={rel[tpv == '24h'].mean():+.3f} vs Ctrl mean={rel[tpv == 'Ctrl'].mean():+.3f}"
    f"  → 对比旧口径 -0.218 vs -0.401, p=0.0734")

# ================= S1: 剂量梯度 JT =================
log("="*76); log(f"S1 | 剂量梯度 naive_A < rCHI < CCI_A < CCI+HS (Jonckheere) [ver={VER}]")
def jt(groups_vals):
    # JT = Σ_{i<j} [#(x_i>x_j) + 0.5·#ties]; 递增趋势 → J<EJ → z<0 → p=Φ(z)
    ns = [len(v) for v in groups_vals]
    J = sum(((a[:,None] > b[None,:]).sum() + 0.5 * (a[:,None] == b[None,:]).sum())
            for i, a in enumerate(groups_vals) for b in groups_vals[i+1:])
    N = sum(ns)
    EJ = (N**2 - sum(n*n for n in ns)) / 4
    VJ = (N**2*(2*N+3) - sum(n*n*(2*n+3) for n in ns)) / 72
    z = (J - EJ) / math.sqrt(VJ)
    from statistics import NormalDist
    return J, z, NormalDist().cdf(z)  # 单侧递增: z 越负越显著
gv = lambda samples, col: pm[(pm["sample"].isin(samples)) & (pm.cell_type == "Oligo")][col].values
for col, nm in [("pct_p21","pct_p21+"), ("senmayo","SenMayo")]:
    seq = [np.asarray(gv(s, col), float) for s in (
        ["RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16"],   # naive 子批A池
        ["RNASEQ09","RNASEQ11","RNASEQ12"],               # rCHI
        ["RNASEQ03","RNASEQ06","RNASEQ15"],               # CCI 子批A侧
        ["RNASEQ01","RNASEQ07","RNASEQ13"])]              # CCI+HS
    J, z, p = jt(seq)
    log(f"    {nm:9s}: J={J:.0f} z={z:+.2f} 单侧p={p:.3g} (渐近; 组均值 "
        + " < ".join(f"{v.mean():.3g}" for v in seq) + ")")

# ================= S2: 同鼠配对 ipsi vs contra =================
log("="*76); log(f"S2 | 同鼠配对 17L/18L/19L vs 17R/18R/19R [ver={VER}]")
ol_cnt = df[df.cell_type == "Oligo"].groupby("sample").size()
excl = [f"RNASEQ{s}{sd}" for s in (17, 18, 19) for sd in "LR"
        if ol_cnt.get(f"RNASEQ{s}{sd}", 0) < MIN_CELLS]
if excl:
    log(f"    (Oligo <{MIN_CELLS} 细胞剔除的配对样本: {', '.join(excl)} — 其配对整对不进检验)")
for col, nm in [("pct_p21","pct_p21+"), ("senmayo","SenMayo")]:
    pairs = []
    for s in (17, 18, 19):
        Li = pm[(pm["sample"] == f"RNASEQ{s}L") & (pm.cell_type == "Oligo")]
        Ri = pm[(pm["sample"] == f"RNASEQ{s}R") & (pm.cell_type == "Oligo")]
        if len(Li) and len(Ri):
            pairs.append((s, Li[col].iloc[0], Ri[col].iloc[0]))
    if len(pairs) < 2:
        log(f"    {nm:9s}: 可用配对 {len(pairs)} <2 → 不可评 (NE)")
        continue
    mice = "/".join(str(p[0]) for p in pairs)
    L_ = [p[1] for p in pairs]; R_ = [p[2] for p in pairs]
    try:
        pw = wilcoxon(L_, R_).pvalue if len(L_) >= 3 else float("nan")
    except ValueError:
        pw = float("nan")  # 全零差异时 wilcoxon 报错, 如实标 nan
    diffpct = 100*np.mean([(l-r)/r if r > 0 else np.nan for l, r in zip(L_, R_)])
    log(f"    {nm:9s} [{mice}]: ipsi={['%.3g'%v for v in L_]} contra={['%.3g'%v for v in R_]} "
        f"配对Wilcoxon p={pw:.3g}  平均升幅={diffpct:+.0f}% (预登记: 对侧差<20%→停; 配对n={len(pairs)})")

# ================= S3: 雌性 =================
log("="*76); log(f"S3 | 雌性 CCI_F vs naive_F [ver={VER}] (n=3v3, exact min p=0.1)")
for col, nm in [("pct_p21","pct_p21+"), ("senmayo","SenMayo")]:
    mwu(g("CCI_F","Oligo",col), g("naive_F","Oligo",col), f"雌性{nm}")

log("="*76)
log(f"[s1_05b] DONE ver={VER}")
KEY = dict(version=VER, p1=p1, p1_or=ors,
           p2_n_types=int(len(p2)), p2_n_sig=int(p2.sig.sum()),
           p2_n_clean=int(len(clean)), p2_n_clean_sig=int(clean.sig.sum()),
           p3_p=float(p3p), p3_mean_24h=float(rel[tpv=="24h"].mean()),
           p3_mean_ctrl=float(rel[tpv=="Ctrl"].mean()))
with open(os.path.join(OUTD, f"readout_keynumbers_{VER}.json"), "w") as f:
    json.dump(KEY, f, indent=1, default=float)
with open(os.path.join(OUTD, f"readout_summary_{VER}.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
