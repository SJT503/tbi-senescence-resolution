# s1_14_gse269748_7d6mo.R — GSE269748 7d/6mo: 清洗 + Cdkn1a/Cdkn2a 鼠级检出率
# 目的: 把急性瞬态从 4 天(GSE277487) 延伸到 6 个月; 回答 Cdkn2a 是否同形状
# 输入:
#   A. 24h/Ctrl 臂: processed_data/s1_24h_pilot/sce_final/  (Step-1 已清洗, floor=1921, 28 样本)
#   B. 7d/6mo 臂:  results/s1_269748_7d6mo/{counts.mtx, genes.txt, cells.csv}  (本次抽取)
# 清洗(7d/6mo 段, 对齐 s1_02→s1_03): decontX → 去双联体(is_dbl) → 共同 floor 下采样
# 统计: 鼠级检出率; 6mo/7d vs 共享 Ctrl —— ⚠️ sham 只在 Ctrl 时点(设计层限制, 必然含 aging 混杂)
# 用法: Rscript s1_14_gse269748_7d6mo.R
# 产物: results/s1_269748_7d6mo/{permouse.csv, LOG.txt, celltable.csv.gz}
suppressPackageStartupMessages({
  library(Matrix); library(decontX); library(SingleCellExperiment); library(Seurat)
})

SEN   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SUBD  <- file.path(SEN, "results/s1_269748_7d6mo")
SCE24 <- file.path(SEN, "processed_data/s1_24h_pilot/sce_final")
OUTD  <- SUBD
SEED  <- 20260925
FLOOR <- 1921L                      # 与 Step-1 一致, 保证跨时点同口径
L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }

# ============ A. 24h / Ctrl 臂 (Step-1 已清洗) ============
log("=== A. 24h/Ctrl 臂 (Step-1 sce_final) ===")
# ⚠️ 样本名含字母后缀 (RNASEQ17L/R, 20L-22L), 禁用 as.integer() 转换 —— 会静默变 NA
sids24 <- sub("\\.rds$", "", list.files(SCE24, pattern = "\\.rds$"))
stopifnot(length(sids24) == 28, !any(is.na(sids24)))
log(sprintf("[A] %d 样本: %s", length(sids24), paste(head(sids24, 6), collapse=", ")))
# A 臂的 sc$group 是分组名, 必须映射到时点 (Ctrl=sham / 24h=injury)
GRP2TP <- c(naive_M="Ctrl", naive_F="Ctrl",
            CCI_ipsi_M="24h", CCI_HS_M="24h", rCHI_M="24h",
            CCI_contra_M="24h", CCI_F="24h")
A <- list()
for (sid in sids24) {
  sc <- readRDS(file.path(SCE24, paste0(sid, ".rds")))
  g <- as.character(sc$group[1])
  A[[length(A)+1]] <- data.frame(
    cell_id = paste0(sid, "_", colnames(sc)),
    sample  = sid,
    time_point = unname(GRP2TP[g]),
    grp_raw = g,
    nUMI = Matrix::colSums(counts(sc)),
    Cdkn1a = as.numeric(counts(sc)["Cdkn1a", ] > 0),
    Cdkn2a = as.numeric(counts(sc)["Cdkn2a", ] > 0),
    stringsAsFactors = FALSE)
}
A <- do.call(rbind, A)
stopifnot(!any(is.na(A$time_point)))
log(sprintf("[A] %d 细胞 | 时点: %s", nrow(A),
            paste(names(table(A$time_point)), table(A$time_point), sep="=", collapse=" ")))

# ============ B. 7d / 6mo 臂 (本次抽取, 需清洗) ============
log("\n=== B. 7d/6mo 臂 (新抽取, 跑 decontX) ===")
M <- as(Matrix::readMM(file.path(SUBD, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(SUBD, "genes.txt")); cells <- read.csv(file.path(SUBD, "cells.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(M) == length(genes), ncol(M) == nrow(cells))
rownames(M) <- genes; colnames(M) <- cells$cell_id
log(sprintf("[B] %d 基因 x %d 细胞", nrow(M), ncol(M)))

# 去双联体 (用 meta 的 is_dbl 标记)
if ("is_dbl" %in% colnames(cells)) {
  keep <- cells$is_dbl == 0
  log(sprintf("[B] 去双联体: 剔除 %d (%.1f%%), 保留 %d",
              sum(!keep), 100*mean(!keep), sum(keep)))
  M <- M[, keep]; cells <- cells[keep, ]
} else log("[B] !! 无 is_dbl 列, 跳过双联体过滤")

# decontX 按样本跑 (与 s1_02 一致)
t0 <- Sys.time()
cleaned <- list()
for (sm in unique(cells$sample)) {
  sel <- which(cells$sample == sm)
  if (length(sel) < 50) { log(sprintf("  %-10s n=%d <50, 跳过", sm, length(sel))); next }
  m0 <- M[, sel, drop = FALSE]
  # ⚠️ decontX 要求 size factors > 0: 必须先滤掉全零细胞 (ScGly_12 曾因此整样本失败)
  nz <- Matrix::colSums(m0) > 0
  if (any(!nz)) {
    log(sprintf("  %-10s 滤除 %d 个全零细胞", sm, sum(!nz)))
    m0 <- m0[, nz, drop = FALSE]
  }
  if (ncol(m0) < 50) { log(sprintf("  %-10s 有效细胞 %d <50, 跳过", sm, ncol(m0))); next }
  m <- m0
  dc <- tryCatch({
    sce <- SingleCellExperiment(list(counts = m))
    sce <- decontX(sce, seed = SEED)
    # ⚠️ decontX 的去污染矩阵在 "decontXcounts" assay; counts(sce) 仍是原始矩阵 (s1_02 正确用法在此)
    list(cnt = round(as(assay(sce, "decontXcounts"), "dgCMatrix")),
         contam = sce$decontX_contamination)
  }, error = function(e) { log(sprintf("  %s decontX 失败: %s", sm, conditionMessage(e))); NULL })
  if (is.null(dc)) next
  cleaned[[sm]] <- data.frame(
    cell_id = colnames(m), sample = sm,
    time_point = cells$time_point[match(colnames(m), cells$cell_id)],
    cell_type  = cells$cell_type[match(colnames(m), cells$cell_id)],
    nUMI = Matrix::colSums(dc$cnt),
    Cdkn1a = as.numeric(dc$cnt["Cdkn1a", ] > 0),
    Cdkn2a = as.numeric(dc$cnt["Cdkn2a", ] > 0),
    contam = dc$contam, stringsAsFactors = FALSE)
  log(sprintf("  %-10s n=%5d nUMI中位=%5.0f 污染中位=%.3f | Cdkn1a %.2f%% Cdkn2a %.2f%%",
              sm, length(sel), median(Matrix::colSums(dc$cnt)), median(dc$contam),
              100*mean(dc$cnt["Cdkn1a", ] > 0), 100*mean(dc$cnt["Cdkn2a", ] > 0)))
}
B <- do.call(rbind, cleaned)
log(sprintf("[B] 清洗后 %d 细胞 (%.1f 分钟)", nrow(B), as.numeric(difftime(Sys.time(), t0, units="mins"))))

# ============ C. 合并 + 共同 floor 下采样 ============
A$cell_type <- NA_character_   # A 需要从 x23 补类型
x23 <- read.csv(gzfile(file.path(SEN, "results/fullrun/x23_scores_GSE269748.csv.gz")),
                stringsAsFactors = FALSE, colClasses = c(cell_id="character"))
A$cell_type <- x23$cell_type[match(A$cell_id, x23$cell_id)]
A$contam <- NA_real_
common <- intersect(colnames(A), colnames(B))
allc <- rbind(A[, common], B[, common])
log(sprintf("\n=== C. 合并 %d 细胞 (A=%d, B=%d)", nrow(allc), nrow(A), nrow(B)))
log("时点分布:")
print(table(allc$time_point, useNA="ifany"))

# 目标细胞型
targets <- c("Microglia","Astrocyte","Oligo","OPC")
allc <- allc[!is.na(allc$cell_type) & allc$cell_type %in% targets, ]
allc <- allc[!is.na(allc$nUMI) & allc$nUMI > 0, ]
log(sprintf("限定 4 型后: %d 细胞", nrow(allc)))

# 共同 floor 下采样 (解析式, 与 s1_12/s1_13 同一函数)
pdet <- function(k, N, F) {
  k <- as.numeric(k); N <- as.numeric(N)
  r <- ifelse(N <= F, as.numeric(k > 0), NA_real_)
  dn <- !(N <= F)
  if (any(dn)) {
    kk <- k[dn]; NN <- N[dn]
    p <- ifelse(kk <= 0, 0, NA_real_)
    nz <- kk > 0
    if (any(nz)) {
      lc <- function(n, r2) lgamma(n+1) - lgamma(r2+1) - lgamma(n-r2+1)
      p[nz] <- 1 - exp(lc(NN[nz]-kk[nz], F) - lc(NN[nz], F))
    }
    r[dn] <- p
  }
  pmin(pmax(r, 0), 1)
}
allc$p21 <- 100 * pdet(allc$Cdkn1a, allc$nUMI, FLOOR)
allc$p16 <- 100 * pdet(allc$Cdkn2a, allc$nUMI, FLOOR)
write.csv(allc, gzfile(file.path(OUTD, "celltable.csv.gz")), row.names = FALSE)

# 鼠级
res <- list()
for (ty in targets) for (g in c("p21","p16")) {
  d <- allc[allc$cell_type == ty, c("sample","time_point", g)]
  ag <- aggregate(d[[g]] ~ sample + time_point, data = d,
                  FUN = function(v) if (length(v) >= 10) mean(v) else NA_real_)
  n  <- aggregate(d[[g]] ~ sample + time_point, data = d, FUN = length)
  ag$n <- n[[3]]; ag$metric <- g; ag$celltype <- ty
  res[[length(res)+1]] <- ag[!is.na(ag[[3]]), ]
}
R <- do.call(rbind, res); names(R)[3] <- "value"
write.csv(R, file.path(OUTD, "permouse.csv"), row.names = FALSE)

# ============ D. 报告 ============
log(paste0("\n", strrep("=", 84)))
log("Cdkn1a (p21) 与 Cdkn2a (p16) 的鼠级检出率 —— 跨时点")
for (ty in targets) {
  log(sprintf("\n########## %s ##########", ty))
  for (g in c("p21","p16")) {
    s <- R[R$celltype == ty & R$metric == g, ]
    if (!nrow(s)) next
    s$time_point <- factor(s$time_point, levels = c("Ctrl","24h","7d","6mo"))
    s <- s[order(s$time_point), ]
    by <- split(s$value, s$time_point)
    log(sprintf("  %s: %s", g, paste(sprintf("%s=%.3f%%(n=%d鼠)", names(by),
                                            sapply(by, mean), sapply(by, length)), collapse=" | ")))
    if ("Ctrl" %in% names(by)) {
      for (tp in setdiff(names(by), "Ctrl")) {
        if (length(by[[tp]]) < 2) next
        p <- tryCatch(wilcox.test(by[[tp]], by$Ctrl, exact=TRUE)$p.value, error=function(e) NA_real_)
        d <- mean(by[[tp]]) - mean(by$Ctrl)
        log(sprintf("     vs Ctrl @ %-4s  差 %+7.3f pp   p=%.4g  (n=%d v %d)", tp, d, p,
                    length(by[[tp]]), length(by$Ctrl)))
      }
    }
  }
}
writeLines(L, file.path(OUTD, "LOG.txt"))
cat("\n[s1_14] DONE\n")
