# s1_15_convergence_probe.R — 判别 7d/6mo 的 "~9% 跨型收敛" 是技术伪影还是生物学
# 背景: s1_14 重跑后发现 7d/6mo 时 Microglia/Astrocyte/Oligo 的 Cdkn1a 检出率全部落在 8.8-9.4%,
#       而 Ctrl 批三型差 14 倍 (14.6/2.65/1.04%). 跨两个文库(ScGly/MB)收敛到同一值 → 可疑
# 三假设:
#   H1 ambient 地板: 中性低表达基因 (Bax/Mcl1) 也跨型收敛; 身份基因保持分叉
#   H2 标签错乱:     身份基因崩塌 (Plp1 in Oligo 从 ~98% 掉到 <60%)
#   H3 生物学:       只有 Cdkn1a 收敛
# A 臂 (Ctrl/24h, Step-1 已清洗): readout_celltable(Cdkn1a/Bax/Mcl1, 已下采样@1921)
#     + _covary_export(Plp1/Mbp/Hexb, 去污染后未下采样 → 解析 pdet@1921)
# B 臂 (7d/6mo): counts.mtx 子集 → 逐样本 decontX(滤全零细胞) → 解析 pdet@1921
# 用法: Rscript s1_15_convergence_probe.R
# 产物: results/s1_269748_7d6mo/{convergence_panel.csv, LOG_CONVERGENCE.txt}
suppressPackageStartupMessages({
  library(Matrix); library(decontX); library(SingleCellExperiment)
})

SEN   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SUBD  <- file.path(SEN, "results/s1_269748_7d6mo")
OUTD  <- SUBD
FLOOR <- 1921L
SEED  <- 20260925
L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }

# ---- 解析下采样检出概率 (与 s1_12/s1_13/s1_14 同一公式, 已对拍验证) ----
pdet <- function(k, N, F) {
  k <- as.numeric(k); N <- as.numeric(N)
  r <- ifelse(N <= F, as.numeric(k > 0), NA_real_)
  dn <- which(!(N <= F))
  if (length(dn)) {
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

GRP2TP <- c(naive_M="Ctrl", naive_F="Ctrl", CCI_ipsi_M="24h", CCI_HS_M="24h",
            rCHI_M="24h", CCI_contra_M="24h", CCI_F="24h")
TYPES <- c("Microglia","Astrocyte","Oligo")

# ============ A 臂 ============
log("=== A. Ctrl/24h 臂 (Step-1 已清洗, 读现有 CSV) ===")
ct  <- read.csv(gzfile(file.path(SEN, "processed_data/s1_24h_pilot/readout_celltable.csv.gz")),
                stringsAsFactors=FALSE, colClasses=c(cell_id="character"))
ct  <- ct[ct$version == "main", c("cell_id","sample","group","p21_umis","bax_umis","mcl1_umis")]
cov <- read.csv(gzfile(file.path(SEN, "processed_data/s1_24h_pilot/_covary_export.csv.gz")),
                stringsAsFactors=FALSE, colClasses=c(cell_id="character"))
# ⚠️ _covary_export 自带 cell_type 列 (s1_08a 导出时就带); 标签对照也从这里建, 不再读外部 labels
stopifnot("cell_type" %in% colnames(cov))
lab <- cov[, c("cell_id","cell_type")]          # 全量标签对照 (含 Ambiguous 等)
cov <- cov[cov$cell_type %in% TYPES, ]
log(sprintf("[A] celltable %d 细胞 | covary(含身份基因) %d 细胞", nrow(ct), nrow(cov)))

# Bax/Mcl1/Cdkn1a: 已物理下采样@1921 → 直接 >0
ct2 <- merge(ct, lab, by="cell_id", all.x=TRUE)
ct2 <- ct2[ct2$cell_type %in% TYPES, ]
ct2$time_point <- unname(GRP2TP[ct2$group])
A_det <- list()
for (g in c("p21_umis","bax_umis","mcl1_umis")) {
  gn <- c(p21_umis="Cdkn1a", bax_umis="Bax", mcl1_umis="Mcl1")[[g]]
  d <- ct2[, c("cell_id","sample","time_point","cell_type")]
  d$det <- as.numeric(ct2[[g]] > 0) * 100
  A_det[[gn]] <- d
}
# 身份基因: covary 去污染后未下采样 → pdet@1921
for (g in c("Plp1","Mbp","Hexb")) {
  d <- cov[, c("cell_id","sample","cell_type")]
  d$sample2 <- sub("_.*$", "", d$cell_id)
  d$grp <- d$sample2
  d$time_point <- unname(GRP2TP[d$grp])
  d$det <- 100 * pdet(cov[[g]], cov$nUMI, FLOOR)
  A_det[[g]] <- d[, c("cell_id","sample2","time_point","cell_type","det")]
  names(A_det[[g]])[2] <- "sample"
}

# ============ B 臂 ============
log("\n=== B. 7d/6mo 臂 (子集 → 逐样本 decontX → pdet@1921) ===")
M <- as(Matrix::readMM(file.path(SUBD, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(SUBD, "genes.txt"))
cells <- read.csv(file.path(SUBD, "cells.csv"), stringsAsFactors=FALSE)
rownames(M) <- genes; colnames(M) <- cells$cell_id
cells <- cells[cells$is_dbl == 0, ]
M <- M[, cells$cell_id]
log(sprintf("[B] 去双联体后 %d 细胞", ncol(M)))

PANEL <- intersect(c("Cdkn1a","Bax","Mcl1","Plp1","Mbp","Hexb","Aqp4","Gfap"), rownames(M))
log(sprintf("[B] 面板基因在场: %s", paste(PANEL, collapse=", ")))

B_rows <- list()
for (sm in unique(cells$sample)) {
  sel <- which(cells$sample == sm)
  m0 <- M[, sel, drop=FALSE]
  nz <- Matrix::colSums(m0) > 0
  if (any(!nz)) { m0 <- m0[, nz, drop=FALSE]; sel <- sel[nz] }
  if (ncol(m0) < 50) { log(sprintf("  %-10s 有效 %d <50, 跳过", sm, ncol(m0))); next }
  dc <- tryCatch({
    sce <- SingleCellExperiment(list(counts = m0))
    sce <- decontX(sce, seed = SEED)
    sce
  }, error = function(e) { log(sprintf("  %s decontX 失败: %s", sm, conditionMessage(e))); NULL })
  if (is.null(dc)) next
  cln <- round(as(assay(dc, "decontXcounts"), "dgCMatrix"))   # ⚠️ 去污染矩阵在 decontXcounts, 非 counts()
  tot <- Matrix::colSums(cln)
  df <- data.frame(cell_id = colnames(m0), sample = sm, stringsAsFactors=FALSE)
  df$time_point <- cells$time_point[match(df$cell_id, cells$cell_id)]
  df$cell_type  <- cells$cell_type[match(df$cell_id, cells$cell_id)]
  for (g in PANEL) df[[g]] <- 100 * pdet(as.numeric(cln[g, ]), tot, FLOOR)
  B_rows[[sm]] <- df
  log(sprintf("  %-10s n=%5d 完成", sm, ncol(m0)))
}
B <- do.call(rbind, B_rows)
B <- B[B$cell_type %in% TYPES, ]
write.csv(B, gzfile(file.path(OUTD, "B_panel_decontX.csv.gz")), row.names=FALSE)  # 先落盘, 汇总段再崩不用重跑
log(sprintf("[B] 合计 %d 细胞 (面板已落盘)", nrow(B)))

# ============ 汇总: 逐时点 x 型 x 基因 ============
log(paste0("\n", strrep("=", 90)))
log("检出率面板 (%) —— 判别收敛性质")
outrows <- list()
add_tp <- function(tp, gene, vals_by_type) {
  for (ty in names(vals_by_type)) {
    outrows[[length(outrows)+1]] <<- data.frame(
      time_point=tp, cell_type=ty, gene=gene,
      mean=mean(vals_by_type[[ty]]), n_mice=length(vals_by_type[[ty]]),
      stringsAsFactors=FALSE)
  }
}
for (gn in c("Cdkn1a","Bax","Mcl1")) {
  d <- A_det[[gn]]
  for (tp in c("Ctrl","24h")) {
    s <- d[d$time_point == tp, ]
    v <- list(); for (ty in TYPES) {
      m <- aggregate(det ~ sample, data=s[s$cell_type==ty,], FUN=function(x) if (length(x)>=10) mean(x) else NA_real_)
      v[[ty]] <- m$det[!is.na(m$det)]
    }
    add_tp(tp, gn, v)
  }
}
for (gn in PANEL) {
  for (tp in c("7d","6mo")) {
    s <- B[B$time_point == tp, ]
    v <- list(); for (ty in TYPES) {
      sb <- s[s$cell_type == ty, c(gn, "sample")]
      names(sb)[1] <- "det"                     # 公式两侧必须同源(子集), 否则长度不等
      ag <- aggregate(det ~ sample, data = sb,
                      FUN = function(x) if (length(x) >= 10) mean(x) else NA_real_)
      v[[ty]] <- ag$det[!is.na(ag$det)]
    }
    add_tp(tp, gn, v)
  }
}
O <- do.call(rbind, outrows)
write.csv(O, file.path(OUTD, "convergence_panel.csv"), row.names=FALSE)

# 打印宽表
for (g in unique(O$gene)) {
  s <- O[O$gene == g, ]
  log(sprintf("\n--- %s ---", g))
  for (tp in c("Ctrl","24h","7d","6mo")) {
    r <- s[s$time_point == tp, ]
    if (!nrow(r)) next
    log(sprintf("  %-4s | %s", tp, paste(
      sprintf("%s=%.2f(n=%d鼠)", r$cell_type, r$mean, r$n_mice), collapse="  ")))
  }
}

# 判别
log(paste0("\n", strrep("=", 90)))
log("判别:")
spread <- function(g, tp) {
  r <- O[O$gene==g & O$time_point==tp, ]
  if (nrow(r) < 2) return(NA_real_)
  max(r$mean) - min(r$mean)
}
for (g in c("Bax","Mcl1")) {
  log(sprintf("  %s 跨型极差: Ctrl=%.1f pp | 7d=%.1f | 6mo=%.1f",
              g, spread(g,"Ctrl"), spread(g,"7d"), spread(g,"6mo")))
}
for (g in c("Plp1","Hexb","Mbp")) {
  r24 <- O[O$gene==g & O$time_point=="24h", ]; r7 <- O[O$gene==g & O$time_point=="7d", ]
  r6 <- O[O$gene==g & O$time_point=="6mo", ]
  pick <- function(r, ty) if (nrow(r[r$cell_type==ty,])) r$mean[r$cell_type==ty][1] else NA
  log(sprintf("  身份 %s: Oligo 24h=%.1f%% 7d=%.1f%% 6mo=%.1f%% | Microglia(若适用) 24h=%.1f%% 7d=%.1f%% 6mo=%.1f%%",
              g, pick(r24,"Oligo"), pick(r7,"Oligo"), pick(r6,"Oligo"),
              pick(r24,"Microglia"), pick(r7,"Microglia"), pick(r6,"Microglia")))
}
writeLines(L, file.path(OUTD, "LOG_CONVERGENCE.txt"))
cat("\n[s1_15] DONE\n")
