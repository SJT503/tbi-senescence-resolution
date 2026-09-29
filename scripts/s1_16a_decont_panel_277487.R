# s1_16a_decont_panel_277487.R — GSE277487 逐样本 decontX + 基因面板导出 (补审计 B1)
# 目的: 之前 s1_10/s1_11b 用原始 counts (decontX 只估分未应用); 本脚本对 12 个样本
#       (各自含全部细胞型 —— ambient 估计需要完整群体) 跑 decontX, 导出去污染后面板
# 面板: p21/p16 + 身份(Plp1/Mbp/Hexb/P2ry12/Tmem119/Aqp4/Gfap/Snap25) + DAM(Trem2/Apoe/Lpl/Cst7/Spp1/Cd9/Clec7a/Itgax)
#       + 慢性程序(Stat1/Irf7) + 凋亡标尺(Bax/Mcl1/Bcl2) + 髓系(Lyz2/C1qa)
# 用法: Rscript s1_16a_decont_panel_277487.R
# 产物: results/s1_277487_decont/panel_decont_allcells.csv.gz + LOG.txt
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(decontX); library(SingleCellExperiment)
})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
CKPT <- file.path(SEN, "results/fullrun/ckpt/ckpt_GSE277487.rds")
OUTD <- file.path(SEN, "results/s1_277487_decont"); dir.create(OUTD, recursive=TRUE, showWarnings=FALSE)
SEED <- 20260925

PANEL <- c("Cdkn1a","Cdkn2a",
           "Plp1","Mbp","Mog","Hexb","P2ry12","Tmem119","Aqp4","Gfap","Snap25",
           "Trem2","Apoe","Lpl","Cst7","Spp1","Cd9","Clec7a","Itgax",
           "Stat1","Irf7","Bax","Mcl1","Bcl2","Lyz2","C1qa")

L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }

t0 <- Sys.time()
log("[s1_16a] 读 ckpt ...")
obj <- readRDS(CKPT)
log(sprintf("[s1_16a] %d 细胞 x %d 基因 (%.1f min)",
            ncol(obj), nrow(obj), as.numeric(difftime(Sys.time(), t0, units="mins"))))
PANEL <- intersect(PANEL, rownames(obj))
miss <- setdiff(c("Cdkn1a","Cdkn2a"), PANEL)
if (length(miss)) stop(paste("核心基因缺失:", paste(miss, collapse=",")))
log(sprintf("[s1_16a] 面板在场 %d 基因", length(PANEL)))

obj$sample <- as.character(obj$sample); obj$tp <- as.character(obj$time_point)
smps <- sort(unique(obj$sample))
stopifnot(length(smps) == 12)
parts <- list()
for (sm in smps) {
  ts <- Sys.time()
  sc <- subset(obj, cells = colnames(obj)[obj$sample == sm])
  cnt <- as(SeuratObject::LayerData(sc[["RNA"]], layer="counts"), "dgCMatrix")
  nz <- Matrix::colSums(cnt) > 0
  if (any(!nz)) { cnt <- cnt[, nz]; log(sprintf("  %-24s 滤 %d 全零细胞", sm, sum(!nz))) }
  sce <- SingleCellExperiment(list(counts = cnt))
  sce <- decontX(sce, seed = SEED)
  cln <- round(as(assay(sce, "decontXcounts"), "dgCMatrix"))     # ⚠️ 去污染矩阵在此 assay (教训 A2)
  cln@x[cln@x < 0] <- 0
  # 断言: assay 替换确实生效 (总 UMI 应低于原始)
  dU <- sum(cln@x); dR <- sum(cnt@x)
  stopifnot(dU <= dR)
  cid <- colnames(cnt)
  df <- data.frame(cell_id = cid, sample = sm,
                   time_point = sc$tp[match(cid, colnames(sc))],
                   cell_type  = as.character(sc$cell_type)[match(cid, colnames(sc))],
                   nUMI_raw = as.numeric(Matrix::colSums(cnt)),
                   nUMI = as.numeric(Matrix::colSums(cln)),
                   contam = as.numeric(sce$decontX_contamination),
                   stringsAsFactors = FALSE)
  pg <- as.matrix(cln[PANEL, ])
  pg <- t(pg); colnames(pg) <- PANEL
  df <- cbind(df, as.data.frame(pg))
  parts[[sm]] <- df
  log(sprintf("  %-24s n=%6d | raw UMI %.2e → clean %.2e (-%.1f%%) | 中位污染 %.3f | Cdkn1a %.2f%%→%.2f%% | %.1f min",
              sm, ncol(cnt), dR, dU, 100*(1-dU/dR), median(sce$decontX_contamination),
              100*mean(cnt["Cdkn1a",] > 0), 100*mean(cln["Cdkn1a",] > 0),
              as.numeric(difftime(Sys.time(), ts, units="mins"))))
  rm(sc, cnt, sce, cln, pg); gc(verbose = FALSE)
}
out <- do.call(rbind, parts)
write.csv(out, gzfile(file.path(OUTD, "panel_decont_allcells.csv.gz")), row.names=FALSE)
log(sprintf("\n[s1_16a] DONE: %d 细胞 x %d 列 → %s (%.1f min 总计)",
            nrow(out), ncol(out), OUTD, as.numeric(difftime(Sys.time(), t0, units="mins"))))
writeLines(L, file.path(OUTD, "LOG.txt"))
