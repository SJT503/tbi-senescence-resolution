# s1_24c — atlas UMAP + MG 状态空间嵌入 (scanpy 稀疏路径; Seurat 稠密版 OOM 后替代)
import scanpy as sc, pandas as pd, numpy as np, scipy.io as sio, scipy.sparse as sp
import os, glob, warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

SEN = "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUT = f"{SEN}/results/s1_atlas_gapfill"
LOG = f"{OUT}/LOG_UMAP.txt"
log = lambda *a: print(*a, file=open(LOG, "a", encoding="utf-8"), flush=True)

# ============ ① 24h 批 atlas (28 样本 ~240k) ============
MD = f"{OUT}/mtx"
if os.path.exists(f"{OUT}/atlas_umap_coords.csv.gz"):
    log("[24c] atlas 已存在, 跳过")
else:
    genes = pd.read_csv(f"{MD}/genes.txt", header=None)[0].tolist()
    labels = pd.read_csv(f"{SEN}/processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")
    adatas = []
    for mtxf in sorted(glob.glob(f"{MD}/*.mtx")):
        sm = os.path.basename(mtxf).replace(".mtx", "")
        X = sp.csr_matrix(sio.mmread(mtxf).T)  # cells x genes
        bc = pd.read_csv(f"{MD}/{sm}_barcodes.txt", header=None)[0].tolist()
        a = sc.AnnData(X.astype(np.float32), obs=pd.DataFrame(index=bc),
                       var=pd.DataFrame(index=genes))
        adatas.append(a)
        log(f"[24c] 读入 {sm}: {a.shape}")
    import anndata as ad, gc
    A = ad.concat(adatas, join="inner")
    del adatas; gc.collect()
    A.obs = labels.set_index("cell_id").loc[A.obs_names]
    log(f"[24c] 合并: {A.shape}")
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A); gc.collect()
    sc.pp.highly_variable_genes(A, n_top_genes=3000, flavor="seurat", batch_key="sample")
    sc.pp.pca(A, n_comps=30, use_highly_variable=True, svd_solver="arpack")
    sc.external.pp.harmony_integrate(A, "sample", basis="X_pca", adjusted_basis="X_pca_harmony")
    sc.pp.neighbors(A, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(A, random_state=20260927)
    co = pd.DataFrame({"cell_id": A.obs_names, "umap_1": A.obsm["X_umap"][:,0],
                       "umap_2": A.obsm["X_umap"][:,1]})
    co = co.merge(labels, on="cell_id")
    co.to_csv(f"{OUT}/atlas_umap_coords.csv.gz", index=False, compression="gzip")
    log(f"[24c] atlas 完成: {len(co)} 细胞 → atlas_umap_coords.csv.gz")

# ============ ② MG 状态空间 (L2 输入复用: 4万 MG, Ctrl/24h/7d/6mo) ============
if os.path.exists(f"{OUT}/mg_statespace_coords.csv.gz"):
    log("[24c] MG 嵌入已存在, 跳过")
else:
    L2 = f"{SEN}/results/s1_L2_origin"
    cells = pd.read_csv(f"{L2}/cells.csv", index_col="cell_id")
    X = sp.csr_matrix(sio.mmread(f"{L2}/counts.mtx").T)
    genes2 = pd.read_csv(f"{L2}/genes.txt", header=None)[0].tolist()
    M = sc.AnnData(X.astype(np.float32), obs=cells, var=pd.DataFrame(index=genes2))
    log(f"[24c] MG: {M.shape}")
    sc.pp.normalize_total(M, target_sum=1e4); sc.pp.log1p(M)
    sc.pp.highly_variable_genes(M, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(M, n_comps=30, svd_solver="arpack")
    sc.external.pp.harmony_integrate(M, "lib", basis="X_pca", adjusted_basis="X_pca_harmony")
    sc.pp.neighbors(M, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(M, random_state=20260925)
    co = pd.DataFrame({"cell_id": M.obs_names, "umap_1": M.obsm["X_umap"][:,0],
                       "umap_2": M.obsm["X_umap"][:,1]})
    co = co.merge(cells.reset_index(), on="cell_id")
    co["p21_pos"] = co.Cdkn1a > 0; co["p16_pos"] = co.Cdkn2a > 0
    co.to_csv(f"{OUT}/mg_statespace_coords.csv.gz", index=False, compression="gzip")
    n21 = int((co.p21_pos & (co.time_point=="24h")).sum()); n16 = int((co.p16_pos & (co.time_point=="6mo")).sum())
    log(f"[24c] MG 嵌入完成: {len(co)} 细胞 | 24h p21+ {n21} | 6mo p16+ {n16} → mg_statespace_coords.csv.gz")
log("[24c] DONE")
print(open(LOG, encoding="utf-8").read()[-1500:])
