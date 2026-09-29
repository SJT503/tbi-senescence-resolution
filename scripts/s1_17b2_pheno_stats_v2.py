# s1_17b2_pheno_stats_v2.py — 表型统计 v2: 逐鼠配对 + 深度分层 (修审计 G1/G2/G5/G7)
# v1 问题: ①池化细胞 Fisher = 伪重复(单元应为鼠) ②未控制细胞深度 ③跨时点对比未标描述性 ④Clec7a 不在矩阵被读成 0% 检出
# 设计: 每鼠每基因检出率 → 配对 Wilcoxon(n鼠) + 方向一致性; 深度敏感性 = 限定共同 nUMI 窗口重跑
# 用法: python s1_17b2_pheno_stats_v2.py
# 产物: results/s1_phenotype/{pheno_v2_results.csv, LOG_PHENO_v2.txt}
import os
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
OUTD = os.path.join(SEN, "results/s1_phenotype")
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

GENES = ["Trem2","Apoe","Lpl","Cst7","Spp1","Cd9","Itgax",
         "P2ry12","Tmem119","Stat1","Irf7","Cdkn1a","Cdkn2a","Bcl2","Mcl1","Bax"]
ABSENT = ["Clec7a"]   # 审计 G7: 不在 GSE269748 基因矩阵, 一律不进表

A = pd.read_csv(os.path.join(OUTD, "mg_panel_24h_v2.csv.gz"))
B = pd.read_csv(os.path.join(OUTD, "mg_panel_6mo_v2.csv.gz"))
log(f"[v2] 24h MG {len(A)} ({A['sample'].nunique()} 鼠) | 6mo MG {len(B)} ({B['sample'].nunique()} 鼠)")

def nUMI_report(d, mask_pos, mask_neg, label):
    log(f"  [{label}] nUMI_total 中位: 阳性 {d.loc[mask_pos,'nUMI_total'].median():.0f} "
        f"vs 阴性 {d.loc[mask_neg,'nUMI_total'].median():.0f} "
        f"(Mann-Whitney 略; 比值 {d.loc[mask_pos,'nUMI_total'].median()/max(d.loc[mask_neg,'nUMI_total'].median(),1):.2f}x)")

def per_mouse_paired(d, pos_mask, neg_mask, label, gene_set=GENES):
    """逐鼠检出率 → 配对 Wilcoxon + 方向一致性 (统计单元=鼠)"""
    rows = []
    mice = sorted(d["sample"].unique())
    for g in gene_set:
        if g in ABSENT: continue
        if g not in d.columns: continue
        pos_v, neg_v = [], []
        for m in mice:
            dm = d[d["sample"] == m]
            p = dm.loc[pos_mask[d["sample"] == m], g]
            n = dm.loc[neg_mask[d["sample"] == m], g]
            if len(p) >= 5 and len(n) >= 20:      # 每鼠最少细胞数
                pos_v.append(100 * (p > 0).mean()); neg_v.append(100 * (n > 0).mean())
        if len(pos_v) < 5:
            rows.append(dict(gene=g, n_mice=len(pos_v), note="鼠数不足"))
            continue
        pos_v, neg_v = np.array(pos_v), np.array(neg_v)
        diff = pos_v - neg_v
        try:
            pw = wilcoxon(pos_v, neg_v, zero_method="wilcox").pvalue if np.any(diff != 0) else 1.0
        except Exception:
            pw = np.nan
        rows.append(dict(gene=g, n_mice=len(pos_v), pos_mean=pos_v.mean(), neg_mean=neg_v.mean(),
                         diff_pp=diff.mean(), direction_consist=f"{(diff>0).sum()}/{len(diff)}",
                         p_wilcoxon=pw))
    r = pd.DataFrame(rows)
    log(f"\n--- {label} (逐鼠配对, n=鼠) ---")
    if len(r):
        log(f"  {'基因':>8s} {'n鼠':>4s} {'阳性%':>8s} {'阴性%':>8s} {'差':>7s} {'方向一致':>8s} {'p(Wilcoxon)':>11s}")
        for _, x in r.iterrows():
            if "note" in x and pd.notna(x.get("note")):
                log(f"  {x.gene:>8s} {int(x.n_mice):4d}  {x.note}"); continue
            log(f"  {x['gene']:>8s} {int(x['n_mice']):4d} {x['pos_mean']:8.1f} {x['neg_mean']:8.1f} {x['diff_pp']:+7.1f} "
                f"{x['direction_consist']:>8s} {x['p_wilcoxon']:11.3g}")
    return r

# ========== 24h: p21+ vs p21− ==========
log("\n" + "="*88)
log("### 24h 小胶质: p21+ vs p21− (逐鼠配对) ###")
pA, nA = A.Cdkn1a > 0, A.Cdkn1a == 0
nUMI_report(A, pA, nA, "24h p21±")
r_full = per_mouse_paired(A, pA, nA, "全深度")

# 深度敏感性: 共同 nUMI 窗口 (全样本 IQR)
q1, q3 = A.nUMI_total.quantile([0.25, 0.75])
log(f"\n  [深度敏感性] 限定 nUMI_total ∈ [{q1:.0f}, {q3:.0f}]")
win = (A.nUMI_total >= q1) & (A.nUMI_total <= q3)
nUMI_report(A, pA & win, nA & win, "24h p21± 窗口内")
r_win = per_mouse_paired(A, (pA & win).values, (nA & win).values, "窗口内")

# ========== 6mo: p16+ vs p16− (n小, 描述性) ==========
log("\n" + "="*88)
log("### 6mo 小胶质: p16+ vs p16− (n小, 描述性 + 深度标注) ###")
pB, nB = B.Cdkn2a > 0, B.Cdkn2a == 0
nUMI_report(B, pB, nB, "6mo p16±")
log(f"  p16+ n={int(pB.sum())} | 逐鼠: " + ", ".join(
    f"{m}:{int(((B['sample']==m)&pB).sum())}" for m in sorted(B['sample'].unique())))
log(f"  {'基因':>8s} {'p16+%':>8s} {'p16−%':>8s} {'差':>7s}")
for g in GENES:
    if g in ABSENT or g not in B.columns: continue
    dp = 100*(B.loc[pB, g] > 0).mean(); dn = 100*(B.loc[nB, g] > 0).mean()
    log(f"  {g:>8s} {dp:8.1f} {dn:8.1f} {dp-dn:+7.1f}")

# ========== 跨时点 (描述性) ==========
log("\n" + "="*88)
log("### 跨时点: 6mo p16+ vs 24h p21+ —— ⚠️ 跨文库+深度未匹配, 仅描述性, 不给 p 值 ###")
log(f"  nUMI 中位: 6mo p16+ {B.loc[pB,'nUMI_total'].median():.0f} vs 24h p21+ {A.loc[pA,'nUMI_total'].median():.0f}")
log(f"  {'基因':>8s} {'6mo p16+%':>10s} {'24h p21+%':>10s}")
for g in GENES:
    if g in ABSENT or g not in B.columns or g not in A.columns: continue
    log(f"  {g:>8s} {100*(B.loc[pB,g]>0).mean():10.1f} {100*(A.loc[pA,g]>0).mean():10.1f}")

r_full.to_csv(os.path.join(OUTD, "pheno_v2_24h_full.csv"), index=False)
r_win.to_csv(os.path.join(OUTD, "pheno_v2_24h_window.csv"), index=False)
with open(os.path.join(OUTD, "LOG_PHENO_v2.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
log(f"\n[s1_17b2] DONE → {OUTD}")
