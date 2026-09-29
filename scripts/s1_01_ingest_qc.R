# s1_01_ingest_qc.R — Step-1 管线第1步: 28样本读入 + QC (预登记: mt<10%, nFeature 下限探查后锁死)
# 用法: Rscript s1_01_ingest_qc.R [smoke]
#   smoke 模式 = 只跑 RNASEQ04(naive) + RNASEQ06(CCI ipsi) 2样本
# 入口门禁(协议§3): 样本数=28 断言 (smoke 模式=2)
# 产物: processed_data/s1_24h_pilot/sce_qc/<sample>.rds  +  qc_summary.csv

suppressPackageStartupMessages({
  library(Matrix); library(SingleCellExperiment)
})

SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
RAW <- file.path(SEN, "raw_data/GSE269748_24h")
OUT <- file.path(SEN, "processed_data/s1_24h_pilot/sce_qc")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- QC 参数 (协议§3: 线粒体<10%; nFeature 下限探查后锁死 → 冒烟后在此处锁值) ----
MT_MAX <- 10        # percent.mt 上限 (%) — 预登记
NFEATURE_MIN <- 200 # 锁死值: snRNA 常规下限, 冒烟分布确认后更新并锁死

# ---- 28 样本分组表 (协议§2, 与 GEO_GSM_map.csv 一致) ----
GROUPS <- list(
  `naive_M` = c("RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25"),
  `CCI_ipsi_M` = c("RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L"),
  `CCI_HS_M` = c("RNASEQ01","RNASEQ07","RNASEQ13"),
  `rCHI_M` = c("RNASEQ09","RNASEQ11","RNASEQ12"),
  `CCI_contra_M` = c("RNASEQ17R","RNASEQ18R","RNASEQ19R"),
  `CCI_F` = c("RNASEQ20L","RNASEQ21L","RNASEQ22L"),
  `naive_F` = c("RNASEQ26","RNASEQ27","RNASEQ28")
)
sample2group <- unlist(lapply(names(GROUPS), function(g) setNames(rep(g, length(GROUPS[[g]])), GROUPS[[g]])))

smoke <- length(commandArgs(trailingOnly = TRUE)) > 0 && commandArgs(trailingOnly = TRUE)[1] == "smoke"
targets <- if (smoke) c("RNASEQ04","RNASEQ06") else names(sample2group)

# 门禁: 全量模式必须 28 样本
stopifnot(smoke || length(targets) == 28)

read10x_manual <- function(gsm, sid) {
  p <- file.path(RAW, sprintf("%s_%s_", gsm, sid))
  mat <- readMM(gzfile(paste0(p, "matrix.mtx.gz")))
  feat <- read.delim(gzfile(paste0(p, "features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  bc <- read.delim(gzfile(paste0(p, "barcodes.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)$V1
  keep <- feat$V3 == "Gene Expression"
  mat <- mat[keep, ]
  gs <- make.unique(feat$V2[keep])
  rownames(mat) <- gs; colnames(mat) <- bc
  as(mat, "dgCMatrix")
}

# GSM 映射 (GEO_GSM_map.csv)
map <- read.csv(file.path(SEN, "results/planA_step47_audit_gapfill/GEO_GSM_map.csv"),
                stringsAsFactors = FALSE)
gsm_of <- setNames(map$gsm, map$sample_id)

summ <- data.frame()
for (sid in targets) {
  gsm <- gsm_of[[sid]]
  cnts <- read10x_manual(gsm, sid)
  ncell_raw <- ncol(cnts)
  mito <- grepl("^mt-", rownames(cnts), ignore.case = TRUE)
  nFeat <- Matrix::colSums(cnts > 0)
  nUMI  <- Matrix::colSums(cnts)
  pct_mt <- if (sum(mito) > 0) Matrix::colSums(cnts[mito, ]) / pmax(nUMI, 1) * 100 else rep(0, ncell_raw)
  keep <- nFeat >= NFEATURE_MIN & pct_mt < MT_MAX
  sce <- SingleCellExperiment(
    assays = list(counts = cnts[, keep]),
    colData = DataFrame(sample = sid, group = sample2group[[sid]],
                        nFeature_raw = nFeat[keep], nUMI_raw = nUMI[keep], pct_mt = pct_mt[keep])
  )
  saveRDS(sce, file.path(OUT, paste0(sid, ".rds")))
  summ <- rbind(summ, data.frame(
    sample = sid, group = sample2group[[sid]], cells_raw = ncell_raw, cells_pass = sum(keep),
    med_nFeat = round(median(nFeat[keep])), med_nUMI = round(median(nUMI[keep])),
    med_pct_mt = round(median(pct_mt[keep]), 2),
    p05_nFeat = round(quantile(nFeat[keep], 0.05)),
    min_nFeat = min(nFeat[keep])))
  cat(sprintf("[s1_01] %-10s raw=%5d pass=%5d (%.1f%%) medFeat=%d medUMI=%d medMT=%.2f\n",
              sid, ncell_raw, sum(keep), 100*sum(keep)/ncell_raw,
              as.integer(median(nFeat[keep])), as.integer(median(nUMI[keep])), median(pct_mt[keep])))
  cat("  nFeat 分位 [1%,5%,25%,50%]:", paste(round(quantile(nFeat, c(.01,.05,.25,.5))), collapse="/"),
      " pct_mt>10% 细胞:", sum(pct_mt >= MT_MAX), "\n")
}
write.csv(summ, file.path(OUT, if (smoke) "qc_summary_smoke.csv" else "qc_summary.csv"), row.names = FALSE)
cat("[s1_01] DONE:", length(targets), "samples ->", OUT, "\n")
