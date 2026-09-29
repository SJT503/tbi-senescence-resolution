# s1_08_myeloid_covary.py — 任务3: p21 阳性 vs myeloid ambient 共变检验
# 问题: Step-1 P1 的 oligo p21 检出率升高 (OR=20.5) 是否为残留髓系 ambient RNA 伪影?
#       已知本数据集源数据广泛存在低水平髓系 ambient (memory: oligo-senmayo-doublet-artifact)
# 判据(预登记式): ①鼠级: %p21+ 与髓系标记检出率 的 Spearman; ②细胞级: 污染分数分层后 p21 是否仍升;
#                 ③去污染后矩阵重算 p21 检出率是否存活 (decontX 已在 s1_02 算过)
# 用法: python s1_08_myeloid_covary.py
# 产物: results/s1_identity_audit/{covary_permouse.csv, covary_celllevel.csv, covary_LOG.txt}
import os, math, gzip
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
SCED = os.path.join(SEN, "processed_data/s1_24h_pilot/sce_clean")
LB   = os.path.join(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")
OUTD = os.path.join(SEN, "results/s1_identity_audit"); os.makedirs(OUTD, exist_ok=True)

L = []
def log(s=""):
    print(s, flush=True); L.append(s)

GROUPS = dict(
    naive_M=["RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25"],
    CCI_ipsi_M=["RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L"],
    CCI_HS_M=["RNASEQ01","RNASEQ07","RNASEQ13"],
    rCHI_M=["RNASEQ09","RNASEQ11","RNASEQ12"],
    CCI_contra_M=["RNASEQ17R","RNASEQ18R","RNASEQ19R"],
    CCI_F=["RNASEQ20L","RNASEQ21L","RNASEQ22L"],
    naive_F=["RNASEQ26","RNASEQ27","RNASEQ28"])
SID2GRP = {s: g for g, v in GROUPS.items() for s in v}
ALL28 = [s for v in GROUPS.values() for s in v]

MYELO = ["C1qa","C1qb","Cx3cr1","P2ry12","Lyz2","Csf1r","Hexb","Ctss"]
OLIGO = ["Plp1","Mbp","Mag","Mog","Cldn11","Mal"]
TARGET = "Cdkn1a"

def cliff(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return (sum((x > b).sum() for x in a) - sum((x < b).sum() for x in b)) / (len(a) * len(b))

def read_rds_matrix(path):
    """读 sce_clean rds → (counts dense-subset dict), 用 R 侧导出更稳; 本函数改为读预导 npz"""
    raise RuntimeError("use export path")

# ---- 用 R 一次性导出所需基因的 counts (避免 python 读 rds) ----
EXPORT = os.path.join(SEN, "processed_data/s1_24h_pilot/_covary_export.csv.gz")
if not os.path.exists(EXPORT):
    log(f"[s1_08] 缺少导出文件 {EXPORT}")
    log("        请先运行 scripts/s1_08a_export_counts.R")
    raise SystemExit(1)

log("[s1_08] 读入导出 counts ...")
df = pd.read_csv(EXPORT)
log(f"  {len(df)} cells, 列: {list(df.columns)[:12]} ...")

# ---- 细胞级: 每个细胞的髓系面板检出 / 污染分数 / p21 ----
df["myeloid_det"] = (df[[c for c in MYELO if c in df.columns]] > 0).any(axis=1)
df["n_myeloid_det"] = (df[[c for c in MYELO if c in df.columns]] > 0).sum(axis=1)
df["n_oligo_det"] = (df[[c for c in OLIGO if c in df.columns]] > 0).sum(axis=1)
df["p21_pos"] = df[TARGET] > 0
df["group"] = df["sample"].map(SID2GRP)

oli = df[df.cell_type == "Oligo"].copy()
log(f"[s1_08] Oligo 细胞 {len(oli)} (14 型中)")

# ================= ① 鼠级共变: %p21+ vs 髓系标记检出率 =================
log("\n" + "=" * 72)
log("① 鼠级共变 (Oligo): %p21+  vs  髓系标记检出率 / 污染分数")
rows = []
for sid, sub in oli.groupby("sample"):
    rows.append(dict(sample=sid, group=SID2GRP.get(sid),
                     n=len(sub),
                     pct_p21=100 * sub.p21_pos.mean(),
                     pct_myeloid_any=100 * sub.myeloid_det.mean(),
                     med_contam=sub.decontX_contamination.median(),
                     med_nUMI=sub.nUMI.median()))
pm = pd.DataFrame(rows).sort_values(["group", "sample"])
pm.to_csv(os.path.join(OUTD, "covary_permouse.csv"), index=False)
log(pm.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

for xcol, lab in [("pct_myeloid_any", "%髓系标记+"), ("med_contam", "中位污染分数"), ("med_nUMI", "中位nUMI")]:
    rho, p = spearmanr(pm.pct_p21, pm[xcol])
    log(f"   Spearman(%p21+, {lab:12s}) = {rho:+.3f}  p={p:.4g}   (n={len(pm)} 鼠)")

# 关键: 只取 24h 两臂 (naive_M + CCI_ipsi_M)
sub24 = pm[pm.group.isin(["naive_M", "CCI_ipsi_M"])]
for xcol, lab in [("pct_myeloid_any", "%髓系标记+"), ("med_contam", "中位污染分数")]:
    rho, p = spearmanr(sub24.pct_p21, sub24[xcol])
    log(f"   [24h两臂] Spearman(%p21+, {lab}) = {rho:+.3f}  p={p:.4g}  (n={len(sub24)} 鼠)")

# 偏相关: 控制污染分数后, 组别 vs %p21+
log("\n   组别效应 (exact MWU):")
x = sub24[sub24.group == "CCI_ipsi_M"].pct_p21.values
y = sub24[sub24.group == "naive_M"].pct_p21.values
if len(x) and len(y):
    p = mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue
    log(f"    原始 %p21+: inj {np.mean(x):.4g} (n={len(x)}) vs sham {np.mean(y):.4g} (n={len(y)})  exact p={p:.4g}  δ={cliff(x,y):+.2f}")

# ================= ② 检出偏倚标尺: 髓系基因自身的鼠级检出率 =================
log("\n" + "=" * 72)
log("② 检出偏倚标尺: 髓系基因在 Oligo 内的鼠级检出率 (若 p21 与之同涨 → 污染嫌疑)")
scale = []
inj_sids = list(sub24.loc[sub24.group == "CCI_ipsi_M", "sample"])
shm_sids = list(sub24.loc[sub24.group == "naive_M", "sample"])
for g in MYELO + [TARGET, "Plp1", "Mbp"]:
    if g not in oli.columns: continue
    r = {}
    for sid, sub in oli.groupby("sample"):
        r[sid] = 100 * (sub[g] > 0).mean()
    scale.append(dict(gene=g,
                      inj=float(np.mean([r[s] for s in inj_sids])),
                      sham=float(np.mean([r[s] for s in shm_sids]))))
sc = pd.DataFrame(scale)
sc["inj_minus_sham"] = sc.inj - sc.sham
log(sc.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

# ================= ③ 污染分数分层后 p21 是否仍升 =================
log("\n" + "=" * 72)
log("③ 污染分数分层: 在低/高污染 Oligo 中分别检验 p21 (24h 两臂)")
q33, q66 = oli.decontX_contamination.quantile([.33, .66])
log(f"   三分位点: q33={q33:.4f}  q66={q66:.4f}")
for name, m in [("低污染(<q33)", oli.decontX_contamination < q33),
                ("中污染(q33-q66)", (oli.decontX_contamination >= q33) & (oli.decontX_contamination < q66)),
                ("高污染(>q66)", oli.decontX_contamination >= q66)]:
    s = oli[m & oli.group.isin(["naive_M", "CCI_ipsi_M"])]
    if not len(s): continue
    a = [100 * sub.p21_pos.mean() for _, sub in s[s.group == "CCI_ipsi_M"].groupby("sample")]
    b = [100 * sub.p21_pos.mean() for _, sub in s[s.group == "naive_M"].groupby("sample")]
    if len(a) >= 1 and len(b) >= 1:
        p = mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue
        log(f"   {name:18s} n=%6d | inj {np.mean(a):7.3f}% (n={len(a)}) vs sham {np.mean(b):7.3f}% (n={len(b)})  p={p:.4g}  δ={cliff(a,b):+.2f}")

# ================= ④ 髓系面板阴性细胞中 p21 检出 =================
log("\n" + "=" * 72)
log("④ 严格子集: 只保留「髓系标记全阴性 且 少突标记阳性」的 Oligo")
strict = oli[(~oli.myeloid_det) & (oli.n_oligo_det >= 2)].copy()
log(f"   严格子集 n={len(strict)} / {len(oli)} ({100*len(strict)/max(len(oli),1):.1f}%)")
s24 = strict[strict.group.isin(["naive_M", "CCI_ipsi_M"])]
a = [100 * sub.p21_pos.mean() for _, sub in s24[s24.group == "CCI_ipsi_M"].groupby("sample")]
b = [100 * sub.p21_pos.mean() for _, sub in s24[s24.group == "naive_M"].groupby("sample")]
log(f"   %p21+: inj {np.mean(a):.4g} (n={len(a)}) vs sham {np.mean(b):.4g} (n={len(b)})")
if len(a) and len(b):
    p = mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue
    log(f"   exact p={p:.4g}  δ={cliff(a,b):+.2f}")
# 池化 Fisher (与 s1_05b G8 同口径)
ia, ib = int(s24[s24.group=="CCI_ipsi_M"].p21_pos.sum()), len(s24[s24.group=="CCI_ipsi_M"])
ca, cb = int(s24[s24.group=="naive_M"].p21_pos.sum()), len(s24[s24.group=="naive_M"])
from scipy.stats import fisher_exact
if ib and cb:
    orv, pf = fisher_exact([[ia, ib-ia], [ca, cb-ca]])
    log(f"   池化 Fisher: {ia}/{ib} ({100*ia/ib:.2f}%) vs {ca}/{cb} ({100*ca/cb:.2f}%)  OR={orv:.2f}  p={pf:.3g}")

out = os.path.join(OUTD, "covary_celllevel.csv")
oli[["sample","group","cell_type","p21_pos","myeloid_det","n_myeloid_det","n_oligo_det",
     "decontX_contamination","nUMI"]].to_csv(out, index=False)
log(f"\n[s1_08] DONE → {OUTD}")
with open(os.path.join(OUTD, "covary_LOG.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
