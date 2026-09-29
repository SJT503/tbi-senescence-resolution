# s1_11_celltype_p21_277487.R — GSE277487 各细胞型 p21 时间轴 (任务A: 6h 峰是否少突特异)
# 目的: 回答"6h p21 峰是少突特有, 还是所有细胞型一起升"
# 口径: 与 s1_10 一致 (共同 floor 下采样 + 鼠级检出率 + KW + 成对 exact MWU)
# 用法: Rscript s1_11_celltype_p21_277487.R
# 产物: results/s1_277487_timecourse/{celltype_p21_permouse.csv, celltype_p21_summary.csv, LOG_CELLTYPE.txt}
suppressPackageStartupMessages({library(Seurat); library(Matrix)})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
CKPT <- file.path(SEN, "results/fullrun/ckpt/ckpt_GSE277487.rds")
OUTD <- file.path(SEN, "results/s1_277487_timecourse"); dir.create(OUTD, recursive=TRUE, showWarnings=FALSE)
SEED <- 20260924
TPS  <- c("Ctrl","group1","group2","group3")
LBL  <- c(Ctrl="Ctrl", group1="6h", group2="2d", group3="4d")

# 够算的型 (每时点 >=3 鼠, 每鼠 >=10 细胞)
TYPES <- c("Oligo","Astrocyte","Microglia","ExN","InN","OPC","Pericyte",
           "Ependymal","Proliferat","Ambiguous_Endothel_ExN","Ambiguous_ChoroidPlx_Pericyte")

L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }
log("[s1_11] 跨细胞型 p21 时间轴 (GSE277487)")

obj <- readRDS(CKPT)
obj$tp <- factor(obj$time_point, levels = TPS)

downsample_sparse <- function(m, target, seed) {
  set.seed(seed)
  tot <- Matrix::colSums(m); keep <- tot > 0
  m <- m[, keep, drop=FALSE]; tot <- tot[keep]
  tgt <- pmin(tot, target)
  ri <- integer(0); ci <- integer(0); xv <- integer(0)
  for (j in seq_len(ncol(m))) {
    idx <- which(m[, j] > 0); if (!length(idx)) next
    v <- rmultinom(1L, size=as.integer(round(tgt[j])), prob=as.numeric(m[idx, j]))[, 1L]
    nz <- v > 0; if (!any(nz)) next
    ri <- c(ri, idx[nz]); ci <- c(ci, rep.int(j, sum(nz))); xv <- c(xv, v[nz])
  }
  Matrix::sparseMatrix(i=ri, j=ci, x=as.numeric(xv), dims=c(nrow(m), ncol(m)),
                       dimnames=list(rownames(m), colnames(m)))
}

permouse <- list(); summary_rows <- list()

for (ty in TYPES) {
  sel <- colnames(obj)[obj$cell_type == ty]
  if (length(sel) < 50) { log(sprintf("  skip %-30s cells=%d", ty, length(sel))); next }
  o <- subset(obj, cells = sel)
  cnt <- as(SeuratObject::LayerData(o[["RNA"]], layer="counts"), "dgCMatrix")
  # 共同 floor = 该型四时点中位的最小值
  med <- tapply(Matrix::colSums(cnt), o$tp, median)
  FLOOR <- as.integer(floor(min(med)))
  cnt_ds <- downsample_sparse(cnt, FLOOR, SEED)
  det <- as.numeric(cnt_ds["Cdkn1a", ] > 0)
  d <- data.frame(mouse = as.character(o$sample), tp = as.character(o$tp),
                  p21 = 100*det, nUMI = Matrix::colSums(cnt_ds), stringsAsFactors=FALSE)
  # 鼠级: 每鼠每时点 >=10 细胞才纳入
  agg <- aggregate(cbind(p21, n) ~ mouse + tp,
                   data = transform(d, n = 1), FUN = function(v) v)
  ag <- do.call(rbind, lapply(split(d, list(d$mouse, d$tp), drop=TRUE), function(x)
    if (nrow(x) >= 10) data.frame(mouse=x$mouse[1], tp=x$tp[1], p21=mean(x$p21), n=nrow(x)) else NULL))
  if (is.null(ag) || !nrow(ag)) { log(sprintf("  skip %-30s (无合格鼠)", ty)); next }
  ag$celltype <- ty; permouse[[ty]] <- ag

  by <- split(ag$p21, ag$tp)
  kw <- tryCatch(kruskal.test(p21 ~ factor(tp), data=ag)$p.value, error=function(e) NA_real_)
  getp <- function(g) tryCatch(wilcox.test(ag$p21[ag$tp==g], ag$p21[ag$tp=="Ctrl"], exact=TRUE)$p.value,
                               error=function(e) NA_real_)
  summary_rows[[ty]] <- data.frame(
    celltype = ty, floor = FLOOR, n_mice_per_tp = paste(sapply(TPS, function(g) sum(ag$tp==g)), collapse="/"),
    Ctrl = mean(by$Ctrl), h6 = mean(by$group1), d2 = mean(by$group2), d4 = mean(by$group3),
    fold_6h = mean(by$group1)/max(mean(by$Ctrl), 1e-9),
    KW_p = kw, p_6h = getp("group1"), p_2d = getp("group2"), p_4d = getp("group3"),
    stringsAsFactors = FALSE)
  log(sprintf("  %-30s floor=%4d | Ctrl %6.3f | 6h %6.3f (%.1fx) | 2d %6.3f | 4d %6.3f | KW p=%.4g",
              ty, FLOOR, mean(by$Ctrl), mean(by$group1),
              mean(by$group1)/max(mean(by$Ctrl),1e-9), mean(by$group2), mean(by$group3), kw))
  rm(o, cnt, cnt_ds); gc(verbose=FALSE)
}
rm(obj); gc(verbose=FALSE)

pm <- do.call(rbind, permouse); sm <- do.call(rbind, summary_rows)
write.csv(pm, file.path(OUTD, "celltype_p21_permouse.csv"), row.names=FALSE)
write.csv(sm, file.path(OUTD, "celltype_p21_summary.csv"), row.names=FALSE)

log("\n" + strrep("=", 78))
log("汇总 (按 6h 倍数降序):")
sm2 <- sm[order(-sm$fold_6h), c("celltype","Ctrl","h6","d2","d4","fold_6h","KW_p")]
print(sm2, row.names=FALSE, digits=4)
log("\n关键判读: Oligo 的 6h 倍数是否显著高于其他型?")
log(sprintf("  Oligo fold_6h = %.2f (rank %d/%d)", sm$fold_6h[sm$celltype=="Oligo"],
            sum(sm$fold_6h > sm$fold_6h[sm$celltype=="Oligo"]) + 1, nrow(sm)))
log(sprintf("  其他型 fold_6h 中位 = %.2f, 范围 %.2f-%.2f",
            median(sm$fold_6h[sm$celltype!="Oligo"]),
            min(sm$fold_6h[sm$celltype!="Oligo"]), max(sm$fold_6h[sm$celltype!="Oligo"])))
writeLines(L, file.path(OUTD, "LOG_CELLTYPE.txt"))
cat("\n[s1_11] DONE\n")
