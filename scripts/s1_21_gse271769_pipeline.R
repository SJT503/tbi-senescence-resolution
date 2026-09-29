# s1_21_gse271769_pipeline.R — GSE271769 解禁管线 (v7 第三验证队列)
# 宪法修订: 2026-09-26 PI 口头批准解禁 (旧禁令基于 planA 星形框架, v7 主角=小胶质 p16)
# 口径: 与 GSE269748 Step-1 同一套 → decontX → scDblFinder → nUMI 共同 floor → 基因级检出
# 产物: results/s1_gse271769/{permouse.csv, LOG.txt, celltable.csv.gz}
# 注意: 本数据集为血管富集制备, 小胶质占 40.7% (56,602 细胞) — 只验小胶质, 不碰少突 (11 个, 不可用)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(decontX); library(SingleCellExperiment)
})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
CKPT <- file.path(SEN, "results/fullrun/ckpt/ckpt_GSE271769.rds")
OUTD <- file.path(SEN, "results/s1_gse271769"); dir.create(OUTD, recursive=TRUE, showWarnings=FALSE)
SEED <- 20260926

L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }

t0 <- Sys.time()
log("[s1_21] 读入 ckpt ...")
obj <- readRDS(CKPT)
cnt_all <- SeuratObject::LayerData(obj[["RNA"]], layer="counts")
mg_mask <- obj$cell_type == "Microglia"
log(sprintf("[s1_21] %d 细胞 | Microglia: %d", ncol(obj), sum(mg_mask)))

# 身份自检 (G3 教训)
hexb <- 100*mean(as.numeric(cnt_all["Hexb", mg_mask] > 0))
snap25 <- 100*mean(as.numeric(cnt_all["Snap25", mg_mask] > 0))
log(sprintf("[s1_21] 身份自检: Hexb+ %.1f%% | Snap25+ %.1f%%", hexb, snap25))
stopifnot(hexb > 80, snap25 < 10)

# ---- 逐样本清洗 (与 Step-1 s1_02 同法) ----
smps <- unique(as.character(obj$sample[mg_mask]))
log(sprintf("\n[s1_21] %d 个小胶质样本", length(smps)))
parts <- list()
for (sm in smps) {
  ts <- Sys.time()
  cells_in_sample <- which(mg_mask & as.character(obj$sample) == sm)
  m <- cnt_all[, cells_in_sample]
  nz <- Matrix::colSums(m) > 0
  if (any(!nz)) m <- m[, nz]
  if (ncol(m) < 50) { log(sprintf("  %-28s n=%d <50 跳过", sm, ncol(m))); next }
  sce <- SingleCellExperiment(list(counts = m))
  sce <- decontX(sce, seed = SEED)
  cln <- round(as(assay(sce, "decontXcounts"), "dgCMatrix"))  # ⚠️ 教训A2: 去污染矩阵在 decontXcounts
  cln@x[cln@x < 0] <- 0
  # scDblFinder (单样本调用, 不传 samples)
  library(scDblFinder)
  sce2 <- scDblFinder(SingleCellExperiment(list(counts = cln)))
  singlet <- colData(sce2)$scDblFinder.class == "singlet"
  cln <- cln[, singlet]
  cid <- colnames(cln)
  tp <- as.character(obj$time_point[cells_in_sample])[match(cid, colnames(m))]
  if (any(is.na(tp))) { log(sprintf("  %s: tp有NA, 排查", sm)); next }
  parts[[sm]] <- data.frame(
    cell_id = cid, sample = sm, time_point = tp[1],
    nUMI = Matrix::colSums(cln),
    Cdkn1a = as.numeric(cln["Cdkn1a", ] > 0),
    Cdkn2a = as.numeric(cln["Cdkn2a", ] > 0),
    stringsAsFactors = FALSE)
  log(sprintf("  %-28s n=%5d singlet | contam=%.3f | Cdkn1a %.1f%% Cdkn2a %.2f%% | %.1f min",
              sm, ncol(cln), median(colData(sce)$decontX_contamination),
              100*mean(cln["Cdkn1a",]>0), 100*mean(cln["Cdkn2a",]>0),
              as.numeric(difftime(Sys.time(), ts, units="mins"))))
  rm(sce, sce2, cln, m); gc(verbose=FALSE)
}

D <- do.call(rbind, parts)
log(sprintf("\n[s1_21] 清洗后 %d 细胞 (%d 样本)", nrow(D), length(unique(D$sample))))

# ---- 深度匹配 (共同 floor, 与 Step-1 同法) ----
med_by_tp <- tapply(D$nUMI, D$time_point, median)
FLOOR <- as.integer(floor(min(med_by_tp)))
log(sprintf("各时点中位 nUMI: %s | 共同 floor = %d",
            paste(names(med_by_tp), as.integer(med_by_tp), sep="=", collapse=" "), FLOOR))

# 解析下采样检出概率 (与 s1_12+ 同一公式)
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
D$p21 <- 100 * pdet(D$Cdkn1a, D$nUMI, FLOOR)
D$p16 <- 100 * pdet(D$Cdkn2a, D$nUMI, FLOOR)
write.csv(D, gzfile(file.path(OUTD, "celltable.csv.gz")), row.names=FALSE)

# ---- 鼠级汇总 + 统计 ----
log(paste0("\n", strrep("=", 76)))
log("GSE271769 小胶质 p21/p16 鼠级检出率 (清洗+深度匹配后)")
res <- list()
for (g in c("p21","p16")) {
  ag <- aggregate(D[[g]] ~ sample + time_point, data=D,
                  FUN=function(v) if(length(v)>=10) mean(v) else NA_real_)
  n <- aggregate(D[[g]] ~ sample + time_point, data=D, FUN=length)
  ag$n <- n[[3]]; ag$metric <- g
  res[[g]] <- ag[!is.na(ag[[3]]), ]
}
R <- do.call(rbind, res); names(R)[3] <- "value"
write.csv(R, file.path(OUTD, "permouse.csv"), row.names=FALSE)

TPS <- c("NoTBI","D1","D7","D28")
for (g in c("p21","p16")) {
  s <- R[R$metric == g, ]
  by <- split(s$value, s$time_point)
  by <- by[intersect(names(by), TPS)]
  log(sprintf("\n%s: %s", g, paste(
    sprintf("%s=%.3f%%(n=%d鼠)", names(by), sapply(by, mean), sapply(by, length)),
    collapse=" | ")))
  if ("NoTBI" %in% names(by)) {
    for (tp in setdiff(names(by), "NoTBI")) {
      x <- by[[tp]]; y <- by$NoTBI
      if (length(x) >= 2 && length(y) >= 2) {
        p <- tryCatch(wilcox.test(x, y, exact=TRUE)$p.value, error=function(e) NA_real_)
        log(sprintf("  vs NoTBI @ %-4s  差 %+7.3f pp   p=%.4g  (n=%d v %d)",
                    tp, mean(x)-mean(y), p, length(x), length(y)))
      }
    }
  }
}
writeLines(L, file.path(OUTD, "LOG.txt"))
cat(sprintf("\n[s1_21] DONE (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
