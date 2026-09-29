# s1_24b — 分样本导出 24h 批清洗计数为 mtx (供 scanpy atlas UMAP; 内存安全: 一次一样本)
suppressPackageStartupMessages({library(SingleCellExperiment)})
SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
ODir <- file.path(SEN, "results/s1_atlas_gapfill/mtx"); dir.create(ODir, showWarnings=FALSE, recursive=TRUE)
LOG <- file.path(SEN, "results/s1_atlas_gapfill/LOG_UMAP.txt")
lab <- read.csv(gzfile(file.path(SEN,"processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")))
cat("[s1_24b] 分样本 mtx 导出开始\n", file=LOG, append=TRUE)
for (sm in sort(unique(lab$sample))) {
  sce <- readRDS(file.path(SEN, sprintf("processed_data/s1_24h_pilot/sce_clean/%s.rds", sm)))
  keep <- paste0(sm, "_", colnames(sce)) %in% lab$cell_id
  m <- counts(sce)[, keep, drop=FALSE]
  colnames(m) <- paste0(sm, "_", colnames(m))
  Matrix::writeMM(m, file=file.path(ODir, paste0(sm, ".mtx")))
  writeLines(colnames(m), file.path(ODir, paste0(sm, "_barcodes.txt")))
  if (sm == sort(unique(lab$sample))[1]) writeLines(rownames(sce), file.path(ODir, "genes.txt"))
  cat(sprintf("[s1_24b] %s: %d 基因 x %d 细胞\n", sm, nrow(m), ncol(m)), file=LOG, append=TRUE)
  rm(sce, m); gc(verbose=FALSE)
}
cat("[s1_24b] 导出完成\n", file=LOG, append=TRUE)
