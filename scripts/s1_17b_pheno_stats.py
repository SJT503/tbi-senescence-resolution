# s1_17b_pheno_stats.py — 表型分析第二步: p21+ vs p16+ 小胶质的状态身份
# 问题: 急性 p21+ MG 与慢性 p16+ MG 是同一程序还是两个状态?
# 输入: results/s1_phenotype/{mg_panel_24h.csv.gz, mg_panel_6mo.csv.gz}
# 模块得分 (mean log1p CPM): DAM / Homeostatic / IFN(Stat1+Irf7) / 髓系身份
# 比较: 24h 批 p21+ vs p21−; 6mo 批 p16+ vs p16−; 6mo p16+ vs 24h p21+
# 用法: python s1_17b_pheno_stats.py
# 产物: results/s1_phenotype/{pheno_cells.csv.gz, LOG_PHENO.txt}
import os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

SEN  = r"DD:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project".replace("DD:", "D:")
OUTD = os.path.join(SEN, "results/s1_phenotype")
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

DAM  = ["Trem2","Apoe","Lpl","Cst7","Spp1","Cd9","Clec7a","Itgax"]
HOME = ["P2ry12","Tmem119"]
IFN  = ["Stat1","Irf7"]
MYID = ["Hexb","Csf1r"] if False else ["Hexb","Lyz2","C1qa"]

def score(d, genes, tag):
    """mean log1p(CPM) per cell for present genes; 深度未匹配 → 用检出率口径另报"""
    present = [g for g in genes if g in d.columns]
    tot = d.nUMI if "nUMI" in d.columns else None
    sub = d[present].astype(float)
    if tot is not None:
        cpm = sub.div(tot, axis=0) * 1e4
        return cpm.apply(np.log1p).mean(axis=1)
    return None

A = pd.read_csv(os.path.join(OUTD, "mg_panel_24h.csv.gz"))
B = pd.read_csv(os.path.join(OUTD, "mg_panel_6mo.csv.gz"))
A["nUMI"] = A[[g for g in A.columns if g not in ("cell_id","sample","arm")]].sum(axis=1)
B["nUMI"] = B[[g for g in B.columns if g not in ("cell_id","sample","arm")]].sum(axis=1)
log(f"[s1_17b] 24h MG {len(A)} | 6mo MG {len(B)}")
log(f"  24h: p21+ {100*(A.Cdkn1a>0).mean():.1f}% | p16+ {100*(A.Cdkn2a>0).mean():.2f}%")
log(f"  6mo: p21+ {100*(B.Cdkn1a>0).mean():.1f}% | p16+ {100*(B.Cdkn2a>0).mean():.2f}% (n={int((B.Cdkn2a>0).sum())} 细胞)")

def compare(d, pos_mask, neg_mask, label, min_pos=8):
    n_pos, n_neg = int(pos_mask.sum()), int(neg_mask.sum())
    log(f"\n--- {label}: +{n_pos} vs −{n_neg} ---")
    if n_pos < min_pos:
        log(f"  !! 阳性细胞 <{min_pos}, 只作描述不作检验")
    det = lambda mask, g: 100*(d.loc[mask, g] > 0).mean()
    rows = []
    for g in DAM + HOME + IFN + ["Cdkn1a","Cdkn2a","Bcl2","Mcl1","Bax"]:
        if g not in d.columns: continue
        dp, dn = det(pos_mask, g), det(neg_mask, g)
        rows.append((g, dp, dn, dp-dn))
    if n_pos >= min_pos:
        log(f"  {'基因':>8s} {'阳性%':>8s} {'阴性%':>8s} {'差':>8s}  Fisher p")
        for g, dp, dn, dd in rows:
            a = int((d.loc[pos_mask, g] > 0).sum()); b = n_pos - a
            c = int((d.loc[neg_mask, g] > 0).sum()); e = n_neg - c
            p = fisher_exact([[a,b],[c,e]][0] if False else [[a,b],[c,e]])[1]
            log(f"  {g:>8s} {dp:8.1f} {dn:8.1f} {dd:+8.1f}  {p:.3g}")
    else:
        log(f"  {'基因':>8s} {'阳性%':>8s} {'阴性%':>8s}")
        for g, dp, dn, dd in rows:
            log(f"  {g:>8s} {dp:8.1f} {dn:8.1f}")
    # 模块检出计数 (每细胞 ≥1 个模块基因即计阳性)
    for name, gs in [("DAM", DAM), ("Homeostatic", HOME), ("IFN", IFN)]:
        present = [g for g in gs if g in d.columns]
        pos_any = (d.loc[pos_mask, present] > 0).any(axis=1)
        neg_any = (d.loc[neg_mask, present] > 0).any(axis=1)
        log(f"  [{name:>11s}] 阳性 {(100*pos_any.mean()):5.1f}% vs 阴性 {(100*neg_any.mean()):5.1f}%")

# 24h: p21+ vs p21−
compare(A, A.Cdkn1a > 0, A.Cdkn1a == 0, "24h 小胶质 p21+ vs p21−")
# 6mo: p16+ vs p16−
compare(B, B.Cdkn2a > 0, B.Cdkn2a == 0, "6mo 小胶质 p16+ vs p16−", min_pos=10)
# 跨时点: 6mo p16+ vs 24h p21+ (同一身份?)
AB = pd.concat([A, B], ignore_index=True)
compare(AB, (AB.arm=="6mo") & (AB.Cdkn2a>0), (AB.arm=="24h") & (AB.Cdkn1a>0),
        "6mo p16+ vs 24h p21+ (跨时点)", min_pos=10)

with open(os.path.join(OUTD, "LOG_PHENO.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
log(f"\n[s1_17b] DONE → {OUTD}")
