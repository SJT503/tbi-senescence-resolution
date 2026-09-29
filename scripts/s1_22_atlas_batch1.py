# s1_22 — 穷尽补缺 Batch 1: 计数分布/阈值敏感性/共表达/组成/pseudobulk (全部从现有 CSV)
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
import os

SEN = "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUT = os.path.join(SEN, "results/s1_atlas_gapfill"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "LOG_BATCH1.txt")
log = lambda *a: print(*a, file=open(LOG, "a", encoding="utf-8"))
log("="*80 + "\n[s1_22] Batch 1: 从现有 CSV 的穷尽补缺\n" + "="*80)

# ---------- 数据载入 ----------
ct69 = pd.read_csv(f"{SEN}/manuscript_package_v7/07_figure_data/GSE269748_7d6mo_celltable.csv.gz")
ct71 = pd.read_csv(f"{SEN}/results/s1_gse271769/celltable.csv.gz")
p77  = pd.read_csv(f"{SEN}/results/s1_277487_decont/panel_decont_allcells.csv.gz")

# ---------- ⑥ 阳性细胞计数分布 ----------
log("\n### ⑥ Cdkn1a/Cdkn2a 阳性细胞计数分布 (k=1/2/3+) ###")
for name, df, genes in [("GSE269748", ct69, ["Cdkn1a","Cdkn2a"]), ("GSE271769(MG)", ct71, ["Cdkn1a","Cdkn2a"])]:
    for g in genes:
        pos = df[df[g] > 0]
        if len(pos) == 0: continue
        dist = pos[g].value_counts().sort_index()
        k1 = dist.get(1,0); k2 = dist.get(2,0); k3 = int(dist[dist.index>=3].sum())
        log(f"  [{name}] {g}: 阳性 n={len(pos)} | k=1: {k1} ({100*k1/len(pos):.1f}%) | k=2: {k2} ({100*k2/len(pos):.1f}%) | k>=3: {k3} ({100*k3/len(pos):.1f}%)")
rows=[]
for name, df in [("GSE269748_4types", ct69), ("GSE271769_MG", ct71)]:
    for g in ["Cdkn1a","Cdkn2a"]:
        for tp in df.time_point.unique():
            pos = df[(df.time_point==tp) & (df[g]>0)]
            rows.append(dict(dataset=name, gene=g, time_point=tp, n_pos=len(pos),
                pct_k1=round(100*(pos[g]==1).mean(),1) if len(pos) else None,
                pct_k2=round(100*(pos[g]==2).mean(),1) if len(pos) else None,
                pct_k3plus=round(100*(pos[g]>=3).mean(),1) if len(pos) else None))
pd.DataFrame(rows).to_csv(f"{OUT}/countdist_positive_cells.csv", index=False)

# ---------- ⑦ 阈值敏感性 (k>=2 / k>=3, 原生深度, 鼠级) ----------
log("\n### ⑦ 阈值敏感性: 检出阈值 k>=1 vs >=2 vs >=3 (原生深度, 逐鼠) ###")
def thresh_permouse(df, gene, kmin, groupcol="time_point", ref="Ctrl"):
    agg = df.assign(det=(df[gene]>=kmin)).groupby([groupcol,"sample"]).det.mean().reset_index()
    out={}
    for tp in agg[groupcol].unique():
        if tp==ref: continue
        x = agg[agg[groupcol]==tp].det.values; y = agg[agg[groupcol]==ref].det.values
        try: p = mannwhitneyu(x, y, alternative="two-sided").pvalue
        except Exception: p=np.nan
        out[tp]=(round(100*x.mean(),3), round(100*y.mean(),3), round(x.mean()/max(y.mean(),1e-9),2), round(p,5))
    return out
rows=[]
for gene, label in [("Cdkn1a","p21"), ("Cdkn2a","p16")]:
    for kmin in [1,2,3]:
        res = thresh_permouse(ct69[ct69.cell_type=="Microglia"], gene, kmin)
        for tp,(m,t,fc,p) in res.items():
            rows.append(dict(dataset="GSE269748_MG", gene=label, kmin=kmin, tp=tp, inj_pct=m, ctrl_pct=t, fold=fc, p=p))
        res2 = thresh_permouse(ct71, gene, kmin, ref="NoTBI")
        for tp,(m,t,fc,p) in res2.items():
            rows.append(dict(dataset="GSE271769_MG", gene=label, kmin=kmin, tp=tp, inj_pct=m, ctrl_pct=t, fold=fc, p=p))
tdf = pd.DataFrame(rows); tdf.to_csv(f"{OUT}/threshold_sensitivity.csv", index=False)
log("  (全表 → threshold_sensitivity.csv; 关键行:)")
for _,r in tdf[(tdf.kmin==2)].iterrows():
    log(f"  {r.dataset} {r.gene} k>=2 @{r.tp}: {r.inj_pct}% vs {r.ctrl_pct}% fold={r.fold} p={r.p}")

# ---------- ⑧ 共表达率 (p21∧p16 双阳) ----------
log("\n### ⑧ p21∧p16 共表达率 (占该型该时点细胞 %) ###")
rows=[]
for name, df, ref in [("GSE269748", ct69, "Ctrl"), ("GSE271769_MG", ct71, "NoTBI")]:
    ct = "Microglia" if "cell_type" in df.columns else None
    for tp in df.time_point.unique():
        sub = df[df.time_point==tp]
        if ct: sub = sub[sub.cell_type=="Microglia"]
        dp = 100*((sub.Cdkn1a>0)&(sub.Cdkn2a>0)).mean()
        p16pos = sub[sub.Cdkn2a>0]
        cd = 100*(p16pos.Cdkn1a>0).mean() if len(p16pos) else np.nan
        rows.append(dict(dataset=name, time_point=tp, n_cells=len(sub),
                         double_pos_pct=round(dp,4), pct_Cdkn1a_pos_among_p16pos=round(cd,1) if cd==cd else None,
                         n_p16pos=len(p16pos)))
co = pd.DataFrame(rows); co.to_csv(f"{OUT}/coexpression_rates.csv", index=False)
log(co.to_string(index=False))

# ---------- ⑨ 组成变化 (277487 全型; CT69 四型) ----------
log("\n### ⑨ 细胞型组成 (每鼠占比) ###")
comp = p77.groupby(["time_point","sample","cell_type"]).size().rename("n").reset_index()
tot = comp.groupby(["time_point","sample"]).n.transform("sum")
comp["frac"] = comp.n/tot
comp.to_csv(f"{OUT}/composition_GSE277487_permouse.csv", index=False)
# KW 检验各型 frac 跨 4 组
from scipy.stats import kruskal
sig=[]
for ctp in comp.cell_type.unique():
    groups=[g.frac.values for _,g in comp[comp.cell_type==ctp].groupby("time_point")]
    if all(len(g)>=3 for g in groups) and len(groups)==4:
        try: p = kruskal(*groups).pvalue
        except Exception: p=np.nan
        means = comp[comp.cell_type==ctp].groupby("time_point").frac.mean()
        if p==p and p<0.05: sig.append((ctp, {k:round(v,4) for k,v in means.items()}, round(p,4)))
log("  GSE277487 组成 KW p<0.05 的型: " + (str(sig) if sig else "无"))
# GSE269748 四型 (CT69, 注意 B 臂仅提取 4 型 → 占比为 4 型内占比)
comp69 = ct69.groupby(["time_point","sample","cell_type"]).size().rename("n").reset_index()
tot69 = comp69.groupby(["time_point","sample"]).n.transform("sum")
comp69["frac"] = comp69.n/tot69
comp69.to_csv(f"{OUT}/composition_GSE269748_4types_permouse.csv", index=False)
for ctp in ["Microglia"]:
    for tp in ["24h","7d","6mo"]:
        x = comp69[(comp69.time_point==tp)&(comp69.cell_type==ctp)].frac.values
        y = comp69[(comp69.time_point=="Ctrl")&(comp69.cell_type==ctp)].frac.values
        try: p = mannwhitneyu(x,y).pvalue
        except Exception: p=np.nan
        log(f"  GSE269748 MG 4型内占比 @{tp}: {x.mean():.3f} vs Ctrl {y.mean():.3f} p={p:.4f}" if p==p else "  NA")

# ---------- ⑩ pseudobulk 佐证 ----------
log("\n### ⑩ pseudobulk: 每鼠平均计数 (Cdkn1a/Cdkn2a mean counts per cell) ###")
rows=[]
for name, df, ref, ctfilter in [("GSE269748_MG", ct69, "Ctrl", "Microglia"),
                                 ("GSE271769_MG", ct71, "NoTBI", None)]:
    sub = df[df.cell_type==ctfilter] if ctfilter else df
    for gene,label in [("Cdkn1a","p21"),("Cdkn2a","p16")]:
        agg = sub.groupby(["time_point","sample"])[gene].mean().reset_index()
        for tp in agg.time_point.unique():
            if tp==ref: continue
            x = agg[agg.time_point==tp][gene].values; y = agg[agg.time_point==ref][gene].values
            try: p = mannwhitneyu(x,y,alternative="two-sided").pvalue
            except Exception: p=np.nan
            rows.append(dict(dataset=name, gene=label, tp=tp, inj_mean=round(x.mean(),6), ctrl_mean=round(y.mean(),6),
                             fold=round(x.mean()/max(y.mean(),1e-12),2), p=round(p,5) if p==p else None))
pb = pd.DataFrame(rows); pb.to_csv(f"{OUT}/pseudobulk_mean_counts.csv", index=False)
log(pb.to_string(index=False))
log("\n[s1_22] Batch1 DONE")
print(open(LOG, encoding="utf-8").read())
