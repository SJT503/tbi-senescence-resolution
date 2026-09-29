# s1_18a_L2_mg_matrix.R — L2 起源裁决第一步: 构建跨时点小胶质矩阵 (清洗+双联体剔除+截断)
# 教训应用(G3): 抽取脚本内重算/对齐索引 + 身份自检; 缺失基因显式记录; 统计单元=鼠
# A 臂: sce_clean (Step-1 已去污染, direct-13: naive_M 7 + ipsi 6)
# B 臂: SUBD 4型子集 → 逐样本 decontX (decontXcounts) → MG (7d 3 + 6mo 4)
# 截断: 每样本 ≤2500 MG (seed 20260925)
# 产物: results/s1_L2_origin/{counts.mtx, genes.txt, cells.csv, LOG.txt}
suppressPackageStartupMessages({
  library(Matrix); library(decontX); library(SingleCellExperiment)
})
SEN   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SCLN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
SUBD  <- file.path(SEN, "results/s1_269748_7d6mo")
OUTD  <- file.path(SEN, "results/s1_L2_origin"); dir.create(OUTD, recursive=TRUE, showWarnings=FALSE)
SEED  <- 20260925; CAP <- 2500L
L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }
GRP2TP <- c(naive_M="Ctrl", CCI_ipsi_M="24h")
DIRECT <- c("RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25",
            "RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L")

lab <- read.csv(gzfile(file.path(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")),
                stringsAsFactors=FALSE, colClasses=c(cell_id="character"))[, c("cell_id","sample","group","cell_type")]
lab <- lab[lab$cell_type == "Microglia" & lab$sample %in% DIRECT, ]
log(sprintf("[A] 标签内 direct-13 MG: %d 细胞", nrow(lab)))

parts <- list(); glist <- list()
set.seed(SEED)
for (sid in DIRECT) {
  f <- file.path(SCLN, paste0(sid, ".rds"))
  if (!file.exists(f)) next
  sce <- readRDS(f)
  cid <- paste0(sid, "_", colnames(sce))
  keep <- cid %in% lab$cell_id
  m <- as(counts(sce), "dgCMatrix")[, keep]          # sce_clean counts = decontXcounts (s1_02)
  cid2 <- cid[keep]
  if (length(cid2) > CAP) { i <- sample(length(cid2), CAP); m <- m[, i]; cid2 <- cid2[i] }
  parts[[sid]] <- list(m = m, cid = cid2, tp = unname(GRP2TP[lab$group[match(cid2[1], lab$cell_id)]]),
                       lib = "RNASEQ")
  glist[[sid]] <- rownames(m)
  log(sprintf("  A %-10s MG %5d (%s)", sid, length(cid2), parts[[sid]]$tp))
  rm(sce, m); gc(verbose=FALSE)
}
genesA <- Reduce(intersect, glist)

# ---- B 臂 ----
M <- as(Matrix::readMM(file.path(SUBD, "counts.mtx")), "CsparseMatrix")
genesB <- readLines(file.path(SUBD, "genes.txt"))
cells <- read.csv(file.path(SUBD, "cells.csv"), stringsAsFactors=FALSE)
rownames(M) <- genesB; colnames(M) <- cells$cell_id
cells <- cells[cells$is_dbl == 0 & cells$cell_type == "Microglia", ]
M <- M[, cells$cell_id]   # ⚠️ 过滤 cells 后必须同步子集 M, 否则 sel 列位错位 (本轮已踩)
log(sprintf("\n[B] 7d/6mo MG (双联体剔除后): %d (M 列同步: %d)", nrow(cells), ncol(M)))
for (sm in sort(unique(cells$sample))) {
  sel <- which(cells$sample == sm)
  m0 <- M[, sel]; nz <- Matrix::colSums(m0) > 0
  if (any(!nz)) { m0 <- m0[, nz]; sel <- sel[nz] }
  sce <- SingleCellExperiment(list(counts = m0))
  sce <- decontX(sce, seed = SEED)
  cln <- round(as(assay(sce, "decontXcounts"), "dgCMatrix")); cln@x[cln@x < 0] <- 0
  cid2 <- colnames(cln)
  if (length(cid2) > CAP) { i <- sample(length(cid2), CAP); cln <- cln[, i]; cid2 <- cid2[i] }
  parts[[sm]] <- list(m = cln, cid = cid2,
                      tp = { tpv <- cells$time_point[match(cid2, cells$cell_id)]
                             stopifnot(!any(is.na(tpv)), length(unique(tpv)) == 1); tpv[1] },
                      lib = ifelse(grepl("^MB", sm), "MB", "ScGly"))
  log(sprintf("  B %-10s MG %5d (%s)", sm, length(cid2), parts[[sm]]$tp))
  rm(sce, cln, m0); gc(verbose=FALSE)
}

genes <- intersect(genesA, genesB)
stopifnot(length(genes) > 20000)
log(sprintf("\n基因交集: %d (A %d / B %d)", length(genes), length(genesA), length(genesB)))

mats <- lapply(parts, function(p) p$m[genes, , drop=FALSE])
BIG <- do.call(cbind, mats)
meta <- do.call(rbind, lapply(names(parts), function(nm) {
  p <- parts[[nm]]
  data.frame(cell_id=p$cid, sample=nm, time_point=p$tp, lib=p$lib, stringsAsFactors=FALSE)
}))
rownames(BIG) <- genes
# 身份自检 (G3 教训): MG 标志物
hexb <- 100*mean(BIG["Hexb",] > 0); p2ry12 <- 100*mean(BIG["P2ry12",] > 0); snap25 <- 100*mean(BIG["Snap25",] > 0)
log(sprintf("身份自检: Hexb+ %.1f%% | P2ry12+ %.1f%% | Snap25+ %.1f%%", hexb, p2ry12, snap25))
stopifnot(hexb > 80, snap25 < 10)
# 状态标注列 (清洗计数, 检出式 + 计数)
for (g in c("Cdkn1a","Cdkn2a","Mcl1","Cst7","Itgax")) meta[[g]] <- as.numeric(BIG[g, ])
Matrix::writeMM(BIG, file.path(OUTD, "counts.mtx"))
writeLines(genes, file.path(OUTD, "genes.txt"))
write.csv(meta, file.path(OUTD, "cells.csv"), row.names=FALSE)
log(sprintf("\n[s1_18a] DONE: %d 基因 x %d 细胞 | 时点: %s", length(genes), ncol(BIG),
            paste(names(table(meta$time_point)), table(meta$time_point), sep="=", collapse=" ")))
log(sprintf("  p21+ 率: %s", paste(sapply(c("Ctrl","24h","7d","6mo"), function(tp)
  sprintf("%s=%.1f%%", tp, 100*mean(meta$Cdkn1a[meta$time_point==tp] > 0))), collapse=" ")))
log(sprintf("  p16+ 率: %s", paste(sapply(c("Ctrl","24h","7d","6mo"), function(tp)
  sprintf("%s=%.2f%%", tp, 100*mean(meta$Cdkn2a[meta$time_point==tp] > 0))), collapse=" ")))
writeLines(L, file.path(OUTD, "LOG.txt"))
