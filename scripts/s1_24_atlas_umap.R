# s1_24 — 24h 批 atlas UMAP (28样本 ~240k细胞, 全细胞型) + marker dot plot 数据
suppressPackageStartupMessages({library(Seurat); library(SingleCellExperiment)})
SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUT <- file.path(SEN, "results/s1_atlas_gapfill"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
LOG <- file.path(OUT, "LOG_UMAP.txt"); log <- function(...) cat(...,"\n",file=LOG,append=TRUE)
cat("",file=LOG); set.seed(20260927)

lab <- read.csv(gzfile(file.path(SEN,"processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")))
log(sprintf("[s1_24] 标签: %d 细胞", nrow(lab)))

# marker panel (dot plot 用)
MARKERS <- c("Hexb","P2ry12","Tmem119","Csf1r","Aqp4","Gfap","Slc1a2","Plp1","Mbp","Mog",
             "Pdgfra","Vcan","Snap25","Syt1","Camk2a","Gad1","Gad2","Cldn5","Pecam1",
             "Pdgfrb","Rgs5","Foxj1","Ttr","Lyz2","C1qa","S100a8","Cd3e","Top2a","Mki67","Cdkn1a","Cdkn2a")

sobjs <- list(); marker_parts <- list()
for (sm in sort(unique(lab$sample))) {
  f <- file.path(SEN, sprintf("processed_data/s1_24h_pilot/sce_clean/%s.rds", sm))
  sce <- readRDS(f)
  m <- as.matrix(counts(sce)[intersect(rownames(sce), MARKERS), , drop=FALSE])
  keep <- paste0(sm, "_", colnames(sce)) %in% lab$cell_id
  # dot plot 数据: 每型每基因 检出率+平均表达(阳性细胞)
  mp <- lab[lab$sample==sm, c("cell_id","cell_type")]
  bc <- paste0(sm, "_", colnames(sce))
  idx <- match(mp$cell_id, bc)
  ok <- !is.na(idx)
  mm <- m[, idx[ok], drop=FALSE]
  q <- data.frame(t(mm) > 0); q$cell_type <- mp$cell_type[ok]
  ag <- aggregate(. ~ cell_type, data=q, FUN=mean)
  marker_parts[[sm]] <- ag
  # UMAP 用全表达
  keepcm <- counts(sce)[, keep, drop=FALSE]
  colnames(keepcm) <- paste0(sm, "_", colnames(keepcm))
  sobjs[[sm]] <- CreateSeuratObject(counts=keepcm, project=sm)
  rm(sce, m, keepcm); gc(verbose=FALSE)
  log(sprintf("[s1_24] %s 读入完成", sm))
}
# marker dot plot 汇总 (跨样本平均)
allm <- do.call(rbind, marker_parts)
dot <- aggregate(. ~ cell_type, data=allm, FUN=mean)
write.csv(dot, file.path(OUT, "marker_dotplot_detection.csv"), row.names=FALSE)
log("[s1_24] marker dot plot 数据完成")

# 合并 + UMAP
mb <- merge(sobjs[[1]], sobjs[-1])
mb[["RNA"]] <- NormalizeData(mb[["RNA"]], verbose=FALSE)
mb <- FindVariableFeatures(mb, nfeatures=3000, verbose=FALSE)
mb <- ScaleData(mb, vars.to.regress=NULL, verbose=FALSE)
mb <- RunPCA(mb, npcs=30, verbose=FALSE)
mb <- RunUMAP(mb, dims=1:30, verbose=FALSE)
log("[s1_24] UMAP 完成")

emb <- Embeddings(mb, "umap")
meta <- lab[match(rownames(emb), lab$cell_id), c("cell_id","sample","group","cell_type","method","confidence")]
out <- data.frame(cell_id=rownames(emb), umap_1=emb[,1], umap_2=emb[,2], meta[,-1])
write.csv(out, gzfile(file.path(OUT,"atlas_umap_coords.csv.gz")), row.names=FALSE)
log(sprintf("[s1_24] 坐标导出: %d 细胞 | 类型数: %d", nrow(out), length(unique(out$cell_type))))
log("[s1_24] DONE")
