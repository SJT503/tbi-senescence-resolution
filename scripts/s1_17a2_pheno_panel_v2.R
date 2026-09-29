# s1_17a2_pheno_panel_v2.R — 表型面板 v2: 补总 UMI + 缺失基因显式记录 (修审计 G1/G4/G7)
# v1 问题: 未导出总 nUMI (深度混杂无法回头校正); 缺失基因静默填 0 (Clec7a 被误读为"0%检出")
# 用法: Rscript s1_17a2_pheno_panel_v2.R
# 产物: results/s1_phenotype/{mg_panel_24h_v2.csv.gz, mg_panel_6mo_v2.csv.gz, LOG_v2.txt}
suppressPackageStartupMessages({
  library(Matrix); library(decontX); library(SingleCellExperiment)
})
SEN   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SCLN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
SUBD  <- file.path(SEN, "results/s1_269748_7d6mo")
OUTD  <- file.path(SEN, "results/s1_phenotype")
SEED  <- 20260925
L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }
PANEL <- c("Cdkn1a","Cdkn2a","Trem2","Apoe","Lpl","Cst7","Spp1","Cd9","Clec7a","Itgax",
           "Stat1","Irf7","Hexb","P2ry12","Tmem119","Lyz2","C1qa","Bcl2","Mcl1","Bax","Plp1")

# ---- Part A: 24h 批 direct-13 ----
log("=== A. 24h 批 direct-13 (sce_clean, 已去污染) ===")
lab <- read.csv(gzfile(file.path(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")),
                stringsAsFactors=FALSE, colClasses=c(cell_id="character"))[, c("cell_id","cell_type")]
DIRECT <- c("RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25",
            "RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L")
lab <- lab[grepl(paste0("^(", paste(DIRECT, collapse="|"), ")_"), lab$cell_id), ]
labMG <- lab$cell_id[lab$cell_type == "Microglia"]
missA <- character(0)
partsA <- list()
for (sid in DIRECT) {
  f <- file.path(SCLN, paste0(sid, ".rds"))
  if (!file.exists(f)) next
  sce <- readRDS(f)
  cid <- paste0(sid, "_", colnames(sce))
  keep <- cid %in% labMG
  if (!any(keep)) next
  cnt <- as(counts(sce), "dgCMatrix")
  m <- cnt[, keep]; cid2 <- cid[keep]
  pn <- intersect(PANEL, rownames(m))
  missA <- union(missA, setdiff(PANEL, pn))
  df <- data.frame(cell_id=cid2, sample=sid, arm="24h",
                   nUMI_total=as.numeric(Matrix::colSums(m)), stringsAsFactors=FALSE)
  for (g in PANEL) df[[g]] <- if (g %in% pn) as.numeric(m[g, ]) else 0
  partsA[[sid]] <- df
  rm(sce, cnt, m); gc(verbose=FALSE)
}
A <- do.call(rbind, partsA)
write.csv(A, gzfile(file.path(OUTD, "mg_panel_24h_v2.csv.gz")), row.names=FALSE)
log(sprintf("[A] 24h MG %d 细胞 | 矩阵缺席基因: %s", nrow(A),
            if (length(missA)) paste(missA, collapse=",") else "无"))

# ---- Part B: 6mo MB (decontX) ----
log("\n=== B. 6mo 批 MB (decontX, decontXcounts) ===")
M <- as(Matrix::readMM(file.path(SUBD, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(SUBD, "genes.txt"))
cells <- read.csv(file.path(SUBD, "cells.csv"), stringsAsFactors=FALSE)
rownames(M) <- genes; colnames(M) <- cells$cell_id
cells <- cells[cells$is_dbl == 0, ]; M <- M[, cells$cell_id]
MB <- sort(unique(cells$sample[grepl("^MB", cells$sample)]))
missB <- character(0)
partsB <- list()
for (sm in MB) {
  sel <- which(cells$sample == sm)
  m0 <- M[, sel]; nz <- Matrix::colSums(m0) > 0
  if (any(!nz)) { m0 <- m0[, nz]; sel <- sel[nz] }
  sce <- SingleCellExperiment(list(counts = m0))
  sce <- decontX(sce, seed = SEED)
  cln <- round(as(assay(sce, "decontXcounts"), "dgCMatrix")); cln@x[cln@x < 0] <- 0
  isMG <- cells$cell_type[match(colnames(m0), cells$cell_id)] == "Microglia"
  mm <- cln[, isMG]; cid2 <- colnames(m0)[isMG]
  pn <- intersect(PANEL, rownames(mm))
  missB <- union(missB, setdiff(PANEL, pn))
  df <- data.frame(cell_id=cid2, sample=sm, arm="6mo",
                   nUMI_total=as.numeric(Matrix::colSums(mm)), stringsAsFactors=FALSE)
  for (g in PANEL) df[[g]] <- if (g %in% pn) as.numeric(mm[g, ]) else 0
  partsB[[sm]] <- df
  log(sprintf("  %-8s MG %d | Cdkn2a+ %d | 中位nUMI %.0f", sm, ncol(mm), sum(mm["Cdkn2a",]>0), median(Matrix::colSums(mm))))
  rm(sce, cln, m0, mm); gc(verbose=FALSE)
}
B <- do.call(rbind, partsB)
write.csv(B, gzfile(file.path(OUTD, "mg_panel_6mo_v2.csv.gz")), row.names=FALSE)
log(sprintf("[B] 6mo MG %d 细胞 | 矩阵缺席基因: %s", nrow(B),
            if (length(missB)) paste(missB, collapse=",") else "无"))
writeLines(L, file.path(OUTD, "LOG_v2.txt"))
cat("\n[s1_17a2] DONE\n")
