# s1_02_decontx_dblfinder.R — Step-1 管线第2步: decontX 去 ambient + scDblFinder 去双联体 (逐样本)
# 协议§3 顺序: decontX(filtered矩阵直接用) → scDblFinder(真实工具替代表达阈值代理)
# 用法: Rscript s1_02_decontx_dblfinder.R [smoke]
# 产物: processed_data/s1_24h_pilot/sce_clean/<sample>.rds (counts=round(decontX), 仅 singlet) + decont_summary.csv
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(celda); library(scDblFinder); library(Matrix)
})

SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
IN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_qc")
OUT <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

smoke <- length(commandArgs(trailingOnly = TRUE)) > 0 && commandArgs(trailingOnly = TRUE)[1] == "smoke"
targets <- if (smoke) c("RNASEQ04","RNASEQ06") else
  sub("\\.rds$", "", list.files(IN, pattern = "\\.rds$"))
stopifnot(smoke || length(targets) == 28)

summ <- data.frame()
for (sid in targets) {
  t0 <- Sys.time()
  sce <- readRDS(file.path(IN, paste0(sid, ".rds")))
  n_in <- ncol(sce)
  # --- decontX (ambient 去污染) ---
  set.seed(20260924)
  sce <- decontX(sce)
  dcx <- round(as(assay(sce, "decontXcounts"), "dgCMatrix"))  # 保持 double: scDblFinder 拒收 integer
  # ambient 污染度诊断: decontX_contamination per cell
  contam <- sce$decontX_contamination
  # --- scDblFinder (在去污染计数上) ---
  sce2 <- SingleCellExperiment(assays = list(counts = dcx),
                               colData = DataFrame(sample = sid, group = sce$group,
                                                   nFeature_raw = sce$nFeature_raw, nUMI_raw = sce$nUMI_raw,
                                                   pct_mt = sce$pct_mt, decontX_contamination = contam))
  set.seed(20260924)
  sce2 <- scDblFinder(sce2)
  singlet <- sce2$scDblFinder.class == "singlet"
  sce2 <- sce2[, singlet]
  saveRDS(sce2, file.path(OUT, paste0(sid, ".rds")))
  n_out <- ncol(sce2)
  summ <- rbind(summ, data.frame(
    sample = sid, cells_in = n_in, cells_singlet = n_out,
    dbl_pct = round(100 * (1 - n_out / n_in), 2),
    med_contam = round(median(contam), 4),
    sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)))
  cat(sprintf("[s1_02] %-10s in=%5d singlet=%5d dbl=%.1f%% medContam=%.3f (%.0fs)\n",
              sid, n_in, n_out, 100 * (1 - n_out / n_in), median(contam),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
write.csv(summ, file.path(OUT, if (smoke) "decont_summary_smoke.csv" else "decont_summary.csv"), row.names = FALSE)
cat("[s1_02] DONE:", length(targets), "samples ->", OUT, "\n")
