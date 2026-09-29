# s1_08a_export_counts.R — 为 s1_08_myeloid_covary.py 导出所需基因的 counts 长表
# 输出: processed_data/s1_24h_pilot/_covary_export.csv.gz
# 列: cell_id, sample, cell_type, decontX_contamination, nUMI, <MYELO>, <OLIGO>, Cdkn1a
# ⚠️ 逐样本追加写 (v2): v1 用 rbind 全量拼装 24 万行 + gzip, 内存峰值被杀 (exit 0 但无产物)
suppressPackageStartupMessages({library(SingleCellExperiment); library(Matrix)})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SCED <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
LB   <- file.path(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")
PARTS <- file.path(SEN, "processed_data/s1_24h_pilot/_covary_parts"); dir.create(PARTS, showWarnings = FALSE)
OUT  <- file.path(SEN, "processed_data/s1_24h_pilot/_covary_export.csv.gz")

GENES <- c("C1qa","C1qb","Cx3cr1","P2ry12","Lyz2","Csf1r","Hexb","Ctss",
           "Plp1","Mbp","Mag","Mog","Cldn11","Mal","Cdkn1a")

lab <- read.csv(gzfile(LB), stringsAsFactors = FALSE)[, c("cell_id","cell_type")]
sids <- sort(unique(sub("_.*$", "", lab$cell_id)))
cat(sprintf("[s1_08a] %d 样本, 锚标签 %d 细胞\n", length(sids), nrow(lab)))

for (sid in sids) {
  pf <- file.path(PARTS, paste0(sid, ".csv.gz"))
  if (file.exists(pf)) { cat(sprintf("  %-10s [skip] 已有\n", sid)); next }
  f <- file.path(SCED, paste0(sid, ".rds"))
  if (!file.exists(f)) { cat("  !! 缺", sid, "\n"); next }
  sce <- readRDS(f)
  cid <- paste0(sid, "_", colnames(sce))
  ct  <- lab$cell_type[match(cid, lab$cell_id)]
  cnt <- as(counts(sce), "dgCMatrix")
  g   <- intersect(GENES, rownames(cnt))
  sub <- as.matrix(t(cnt[g, , drop = FALSE]))
  df  <- data.frame(
    cell_id = cid, sample = sid, cell_type = ifelse(is.na(ct), "UNMATCHED", ct),
    decontX_contamination = as.numeric(colData(sce)$decontX_contamination),
    nUMI = Matrix::colSums(cnt), stringsAsFactors = FALSE)
  df <- cbind(df, as.data.frame(sub, stringsAsFactors = FALSE))
  for (gg in setdiff(GENES, g)) df[[gg]] <- 0L
  df <- df[, c("cell_id","sample","cell_type","decontX_contamination","nUMI", GENES)]
  write.csv(df, gzfile(pf), row.names = FALSE)
  cat(sprintf("  %-10s cells=%6d → part written\n", sid, nrow(df)))
  rm(sce, cnt, sub, df); gc(verbose = FALSE)
}

# ---- 合并: 流式拼接 (逐文件读写, 不整体驻留) ----
outcon <- gzfile(OUT, "w")
first <- TRUE
for (sid in sids) {
  pf <- file.path(PARTS, paste0(sid, ".csv.gz"))
  if (!file.exists(pf)) next
  con <- gzfile(pf, "r"); ln <- readLines(con, warn = FALSE); close(con)
  if (length(ln) == 0) next
  if (!first) ln <- ln[-1]           # 去掉表头
  writeLines(ln, outcon)
  first <- FALSE
}
close(outcon)
cat(sprintf("[s1_08a] DONE → %s\n", OUT))
