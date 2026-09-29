# s1_17a_pheno_panel.R — 表型分析第一步: 导出小胶质 p21+/p16+ 表型面板
# 问题: 急性 p21+ 小胶质 vs 慢性 p16+ 小胶质是同一群/同一程序, 还是两个不同状态?
# Part A (24h): sce_clean 的 direct-13 样本 (已去污染), 小胶质 p21+ vs p21− 表型
# Part B (6mo): MB 批 4 样本 4型子集 → decontX → 小胶质 p16+ vs p16− 表型
# 用法: Rscript s1_17a_pheno_panel.R
# 产物: results/s1_phenotype/{mg_panel_24h.csv.gz, mg_panel_6mo.csv.gz, LOG.txt}
suppressPackageStartupMessages({
  library(Matrix); library(decontX); library(SingleCellExperiment)
})

SEN   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SCLN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
SUBD  <- file.path(SEN, "results/s1_269748_7d6mo")
OUTD  <- file.path(SEN, "results/s1_phenotype"); dir.create(OUTD, recursive=TRUE, showWarnings=FALSE)
SEED  <- 20260925
L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }

PANEL <- c("Cdkn1a","Cdkn2a","Trem2","Apoe","Lpl","Cst7","Spp1","Cd9","Clec7a","Itgax",
           "Stat1","Irf7","Hexb","P2ry12","Tmem119","Lyz2","C1qa","Bcl2","Mcl1","Bax","Plp1")

# ============ Part A: 24h 批 (sce_clean 已去污染) ============
log("=== A. 24h 批 direct-13 (sce_clean) ===")
lab <- read.csv(gzfile(file.path(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")),
                stringsAsFactors=FALSE, colClasses=c(cell_id="character"))[, c("cell_id","cell_type")]
# 只取 direct-13 (x23 join 命中集 —— 与 s1_14 A 臂口径一致)
DIRECT <- c("RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25",
            "RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L")
lab <- lab[lab$cell_id %in% paste0(rep(DIRECT, 1), "_") | grepl(paste0("^(", paste(DIRECT, collapse="|"), ")_"), lab$cell_id), ]
labMG <- lab$cell_id[lab$cell_type == "Microglia"]
partsA <- list()
for (sid in DIRECT) {
  f <- file.path(SCLN, paste0(sid, ".rds"))
  if (!file.exists(f)) next
  sce <- readRDS(f)
  cid <- paste0(sid, "_", colnames(sce))
  keep <- cid %in% labMG
  if (!any(keep)) next
  cnt <- as(counts(sce), "dgCMatrix")   # sce_clean 的 counts 已是 decontXcounts (s1_02 正确存法)
  m <- cnt[, keep]; cid2 <- cid[keep]
  pn <- intersect(PANEL, rownames(m))
  df <- data.frame(cell_id=cid2, sample=sid, arm="24h", stringsAsFactors=FALSE)
  for (g in pn) df[[g]] <- as.numeric(m[g, ])
  miss <- setdiff(PANEL, pn); for (g in miss) df[[g]] <- 0
  partsA[[sid]] <- df[, c("cell_id","sample","arm", PANEL)]
  rm(sce, cnt, m); gc(verbose=FALSE)
}
A <- do.call(rbind, partsA)
write.csv(A, gzfile(file.path(OUTD, "mg_panel_24h.csv.gz")), row.names=FALSE)
log(sprintf("[A] 24h 小胶质 %d 细胞 (13 direct 样本)", nrow(A)))

# ============ Part B: 6mo 批 (MB, 现跑 decontX) ============
log("\n=== B. 6mo 批 MB (decontX on 4型子集) ===")
M <- as(Matrix::readMM(file.path(SUBD, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(SUBD, "genes.txt"))
cells <- read.csv(file.path(SUBD, "cells.csv"), stringsAsFactors=FALSE)
rownames(M) <- genes; colnames(M) <- cells$cell_id
cells <- cells[cells$is_dbl == 0, ]
M <- M[, cells$cell_id]
MB <- cells$sample[grepl("^MB", cells$sample)] |> unique() |> sort()
partsB <- list()
for (sm in MB) {
  sel <- which(cells$sample == sm)
  m0 <- M[, sel]; nz <- Matrix::colSums(m0) > 0
  if (any(!nz)) { m0 <- m0[, nz]; sel <- sel[nz] }
  sce <- SingleCellExperiment(list(counts = m0))
  sce <- decontX(sce, seed = SEED)
  cln <- round(as(assay(sce, "decontXcounts"), "dgCMatrix"))
  cln@x[cln@x < 0] <- 0
  isMG <- cells$cell_type[match(colnames(m0), cells$cell_id)] == "Microglia"
  mm <- cln[, isMG]; cid2 <- colnames(m0)[isMG]
  pn <- intersect(PANEL, rownames(mm))
  df <- data.frame(cell_id=cid2, sample=sm, arm="6mo", stringsAsFactors=FALSE)
  for (g in pn) df[[g]] <- as.numeric(mm[g, ])
  miss <- setdiff(PANEL, pn); for (g in miss) df[[g]] <- 0
  partsB[[sm]] <- df[, c("cell_id","sample","arm", PANEL)]
  log(sprintf("  %-8s MG %d 细胞 | Cdkn2a+ %d", sm, ncol(mm), sum(mm["Cdkn2a",] > 0)))
  rm(sce, cln, m0, mm); gc(verbose=FALSE)
}
B <- do.call(rbind, partsB)
write.csv(B, gzfile(file.path(OUTD, "mg_panel_6mo.csv.gz")), row.names=FALSE)
log(sprintf("[B] 6mo 小胶质 %d 细胞 (%d MB 样本)", nrow(B), length(MB)))
writeLines(L, file.path(OUTD, "LOG.txt"))
cat("\n[s1_17a] DONE\n")
