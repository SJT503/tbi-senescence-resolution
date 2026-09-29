# -*- coding: utf-8 -*-
"""s1_19_L3_spatial.py — L3 空间层: 慢性微环境是否白质限定 (GSE319409, 7d Visium)
v7 框架下重新设计: p16+ 太稀疏不可直接检测 → 改测微环境
读出:
  ① 小胶质模块评分: WM vs GM × TBI vs Sham (小胶质是否在伤后 WM 富集)
  ② 基因级空间检出: Cdkn1a(p21 波应已退) / Cdkn2a(p16 任何残余) / Mcl1 / Cst7 / Itgax
  ③ 维持信号: Cxcl12 / Csf1 / Mif / Tnf spot 表达
  ④ 小胶质-WM 共定位: 小胶质评分与 oligo 评分的相关性 (TBI vs Sham)
数据: raw_data/GSE319409/ (12 样本 TBI_Veh vs Sham_Veh, 6v6)
产物: results/s1_L3_spatial/{per_spot_all.csv, per_sample_summary.csv, LOG.txt}
"""
import os, glob, time
import numpy as np, pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import wilcoxon, spearmanr

SEN  = r"D:\2026-04-24\new-chat-2\TBI_astrocyte_scst_project\senescence_project"
RAWD = os.path.join(SEN, "raw_data/GSE319409")
OUTD = os.path.join(SEN, "results/s1_L3_spatial"); os.makedirs(OUTD, exist_ok=True)
SEED = 20260925
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

t0 = time.time()

# ---- 模块定义 ----
MODULES = {
    "microglia":  ["Hexb", "Csf1r", "P2ry12", "Tmem119", "Cx3cr1"],
    "oligo_wm":   ["Plp1", "Mbp", "Mog", "Mag", "Cldn11"],
    "dam":        ["Trem2", "Apoe", "Lpl", "Cst7", "Spp1", "Cd9"],
    "ifn":        ["Stat1", "Irf7", "Ifit1", "Isg15"],
}
# 基因级单基因检出
SINGLE_GENES = ["Cdkn1a", "Cdkn2a", "Mcl1", "Cst7", "Itgax",
                "Cxcl12", "Csf1", "Mif", "Tnf", "Adam17", "Osm", "Il1b"]
# 维持信号配体 (来自 L1 NicheNet)

# ---- 读入样本 ----
sample_dirs = sorted(glob.glob(os.path.join(RAWD, "GSM*")))
log(f"[L3] 找到 {len(sample_dirs)} 个 GSM 目录")

# 读 samples.csv (step25 用的分组信息)
meta_f = os.path.join(RAWD, "samples.csv")
if os.path.exists(meta_f):
    smeta = pd.read_csv(meta_f)
    log(f"[L3] samples.csv: {len(smeta)} 行 | 列: {list(smeta.columns)[:6]}")
    # 过滤 TBI_Veh + Sham_Veh
    grp_col = [c for c in smeta.columns if "group" in c.lower() or "condition" in c.lower() or "sample" in c.lower()]
    if grp_col:
        log(f"  分组列: {grp_col}")
        for c in grp_col[:2]:
            log(f"  {c} 值: {smeta[c].unique()}")
else:
    log("[L3] ⚠️ 无 samples.csv — 从 step25 的 spot_scores_all.csv 取分组")

# 优先: 直接用 step25 的 spot_scores (已有完整基础设施)
step25_f = os.path.join(SEN, "results/planA_step25_spatial_gse319409/spot_scores_all.csv")
if os.path.exists(step25_f):
    log("\n[L3] 使用 step25 的 spot_scores_all.csv (已有基础设施)")
    spots = pd.read_csv(step25_f)
    log(f"[L3] {len(spots)} spots | 列: {list(spots.columns)[:10]}...")
    # 但这个文件可能没有我们要的单基因 — 需要检查
    has_cdkn1a = "Cdkn1a" in spots.columns
    has_cdkn2a = "Cdkn2a" in spots.columns
    log(f"  Cdkn1a 在列: {has_cdkn1a} | Cdkn2a 在列: {has_cdkn2a}")
    if not (has_cdkn1a and has_cdkn2a):
        log("  → 需要从原始矩阵重新计算基因级检出")
        spots = None

# ---- 如果需要, 从原始矩阵计算 ----
if spots is None:
    log("\n[L3] 从原始 Visium 矩阵计算...")
    all_spots = []
    for gsm_dir in sample_dirs:
        gsm = os.path.basename(gsm_dir)
        # 读三件套
        mtx_f = glob.glob(os.path.join(gsm_dir, "*matrix.mtx*"))
        feat_f = glob.glob(os.path.join(gsm_dir, "*features.tsv*"))
        bar_f = glob.glob(os.path.join(gsm_dir, "*barcodes.tsv*"))
        coord_f = glob.glob(os.path.join(gsm_dir, "*tissue_positions*"))
        if not (mtx_f and feat_f and bar_f):
            log(f"  {gsm}: 缺文件, 跳过"); continue
        # 读入
        import gzip as gz
        from scipy.io import mmread
        # 处理 .gz 压缩文件
        def _read_mtx(f):
            if f.endswith(".gz"):
                import shutil, tempfile
                with gz.open(f, "rb") as fin, tempfile.NamedTemporaryFile(delete=False, suffix=".mtx") as fout:
                    shutil.copyfileobj(fin, fout); return fout.name
            return f
        mtx_path = _read_mtx(mtx_f[0])
        X = mmread(mtx_path).T.tocsr()
        def _read_tsv(f):
            if f.endswith(".gz"):
                return pd.read_csv(f, sep="\t", header=None, compression="gzip")
            return pd.read_csv(f, sep="\t", header=None)
        feats = _read_tsv(feat_f[0])
        genes_raw = feats[1].values if len(feats.columns) > 1 else feats[0].values
        # 去重基因名 (Visium 常见重复, 保留首个)
        _, uniq_idx = np.unique(genes_raw, return_index=True)
        uniq_idx = np.sort(uniq_idx)
        genes = genes_raw[uniq_idx]
        X = X[:, uniq_idx]
        def _read_lines(f):
            if f.endswith(".gz"):
                return [l.strip() for l in gz.open(f, "rt")]
            return [l.strip() for l in open(f)]
        barcodes = _read_lines(bar_f[0])
        adata = ad.AnnData(X=X.astype(np.float32),
                           obs=pd.DataFrame(index=barcodes),
                           var=pd.DataFrame(index=pd.Index(genes)))
        # 归一化
        sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
        # 模块评分 (mean expression)
        for mod_name, mod_genes in MODULES.items():
            present = [g for g in mod_genes if g in adata.var_names]
            if present:
                adata.obs[f"mod_{mod_name}"] = np.asarray(adata[:, present].X.mean(axis=1)).ravel()
            else:
                adata.obs[f"mod_{mod_name}"] = 0.0
                log(f"  ⚠️ {gsm} 模块 {mod_name} 0 基因在场")
        # 单基因检出 (raw >0)
        # 重新读原始矩阵(未归一化)来判检出
        X_raw = mmread(mtx_f[0]).T.tocsr()
        for g in SINGLE_GENES:
            if g in adata.var_names:
                gi = list(adata.var_names).index(g)
                adata.obs[f"det_{g}"] = (np.asarray(X_raw[:, gi].todense()).ravel() > 0).astype(int)
            else:
                adata.obs[f"det_{g}"] = 0
        adata.obs["gsm"] = gsm
        # 空间坐标 (如果有)
        if coord_f:
            try:
                coord = pd.read_csv(coord_f[0], index_col=0)
                if "pxl_row_in_fullres" in coord.columns:
                    adata.obs["row"] = coord.loc[barcodes, "pxl_row_in_fullres"].values
                    adata.obs["col"] = coord.loc[barcodes, "pxl_col_in_fullres"].values
            except: pass
        df = adata.obs.copy()
        all_spots.append(df)
        log(f"  {gsm}: {len(df)} spots")
    spots = pd.concat(all_spots, ignore_index=True)

# ---- 分组标注 ----
# 从 step25 per_sample_summary 取分组
ps25 = pd.read_csv(os.path.join(SEN, "results/planA_step25_spatial_gse319409/per_sample_summary.csv"))
gsm2grp = dict(zip(ps25["sample"], ps25["group"]))
spots["group"] = spots["gsm"].map(gsm2grp) if "gsm" in spots.columns else spots.get("group")
log(f"\n[L3] 分组: {spots['group'].value_counts().to_dict()}")

# ---- 白质分区 (oligo 上 1/3, 与 step25 口径一致) ----
oligo_col = [c for c in spots.columns if "oligo" in c.lower() and "mod" in c]
if not oligo_col:
    oligo_col = [c for c in spots.columns if "oligo" in c.lower()]
log(f"[L3] oligo 评分列: {oligo_col[:3]}")
if oligo_col:
    oc = oligo_col[0]
    spots["is_wm"] = spots.groupby("gsm")[oc].transform(lambda x: x >= x.quantile(2/3))
    log(f"[L3] WM spots: {spots['is_wm'].sum()}/{len(spots)} ({100*spots['is_wm'].mean():.0f}%)")

# ---- 统计 ----
log("\n" + "=" * 84)
log("L3 空间层结果 (7d, GSE319409)")

# 1) 小胶质在 WM vs GM, TBI vs Sham
log("\n### ① 小胶质模块: WM vs GM × TBI vs Sham ###")
mg_col = [c for c in spots.columns if "microglia" in c and "mod" in c]
if mg_col:
    mc = mg_col[0]
    for grp in ["Sham", "TBI"]:
        for wm in [True, False]:
            s = spots[(spots.group == grp) & (spots.is_wm == wm)]
            gsm_means = s.groupby("gsm")[mc].mean()
            log(f"  {grp:>4s} {'WM' if wm else 'GM'}: mean={gsm_means.mean():.3f} "
                f"(n={len(gsm_means)} 鼠, 逐鼠: {['%.3f' % v for v in gsm_means.values]})")
    # 检验: TBI 的 WM 小胶质 > GM?
    tbi_wm = spots[(spots.group == "TBI") & (spots.is_wm)].groupby("gsm")[mc].mean()
    tbi_gm = spots[(spots.group == "TBI") & (~spots.is_wm)].groupby("gsm")[mc].mean()
    if len(tbi_wm) >= 3:
        from scipy.stats import wilcoxon
        p = wilcoxon(tbi_wm, tbi_gm).pvalue if len(tbi_wm) >= 3 else np.nan
        log(f"  TBI 内 WM>GM Wilcoxon p = {p:.4f} (n={len(tbi_wm)} 对)")

# 2) 基因级空间检出
log("\n### ② 基因级检出 (spot 级 %) ###")
for g in ["Cdkn1a", "Cdkn2a", "Mcl1", "Cst7", "Itgax", "Cxcl12", "Csf1", "Mif"]:
    dc = f"det_{g}"
    if dc not in spots.columns:
        # 如果来自 step25, 可能用不同列名
        alt = [c for c in spots.columns if g.lower() in c.lower()]
        if alt: dc = alt[0]
        else: continue
    for grp in ["Sham", "TBI"]:
        s = spots[spots.group == grp]
        det_rate = 100 * s[dc].mean() if s[dc].dtype in [int, float, bool] else np.nan
        if np.isnan(det_rate):
            # 可能是连续值
            det_rate = 100 * (s[dc] > 0).mean()
        # WM / GM 分层
        wm_det = 100 * (s[s.is_wm][dc] > 0).mean() if s.is_wm.sum() > 0 else np.nan
        gm_det = 100 * (s[~s.is_wm][dc] > 0).mean() if (~s.is_wm).sum() > 0 else np.nan
        log(f"  {g:>8s} {grp:>4s}: 全spot {det_rate:5.1f}% | WM {wm_det:5.1f}% | GM {gm_det:5.1f}%")

# 3) 小胶质-WM 共定位
log("\n### ③ 小胶质-oligo 评分相关性 (TBI vs Sham) ###")
if mg_col and oligo_col:
    for grp in ["Sham", "TBI"]:
        cors = []
        for g, s in spots[spots.group == grp].groupby("gsm"):
            r, _ = spearmanr(s[mc], s[oc])
            cors.append(r)
        log(f"  {grp:>4s}: Spearman(microglia, oligo) mean={np.mean(cors):.3f} "
            f"(n={len(cors)} 鼠, range {np.min(cors):.3f}~{np.max(cors):.3f})")

# 4) 维持信号在 WM 的富集
log("\n### ④ 维持信号配体在 WM vs GM ###")
for g in ["Cxcl12", "Csf1", "Mif", "Tnf"]:
    dc = f"det_{g}"
    if dc not in spots.columns: continue
    for grp in ["Sham", "TBI"]:
        s = spots[spots.group == grp]
        wm_d = 100 * (s[s.is_wm][dc] > 0).mean() if s.is_wm.sum() > 0 else np.nan
        gm_d = 100 * (s[~s.is_wm][dc] > 0).mean() if (~s.is_wm).sum() > 0 else np.nan
        log(f"  {g:>8s} {grp:>4s}: WM {wm_d:5.1f}% vs GM {gm_d:5.1f}% (Δ={wm_d-gm_d:+.1f}pp)")

# 保存
spots.to_csv(os.path.join(OUTD, "per_spot_all.csv"), index=False)
log(f"\n[L3] DONE ({time.time()-t0:.0f}s) → {OUTD}")
with open(os.path.join(OUTD, "LOG.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
