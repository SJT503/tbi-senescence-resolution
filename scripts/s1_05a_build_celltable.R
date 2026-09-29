# s1_05a_build_celltable.R — Step-1 第5步a: 逐细胞指标表 (SenMayo_real 评分 + p21/Bax/Mcl1 检出)
# 评分口径逐位复刻 x23_senescence_score.R:160-168 (方法1 平均标准化表达):
#   lognorm = log1p(count/libsize*1e4); SenMayo_real = colMeans(lognorm[intersect(118基因)])
# 三版本并行 (协议§4 敏感性±20%): main(=s1_03 已存 sce_final) / lo20 / hi20 (从 sce_clean 以固定种子重下采样)
# 产物: processed_data/s1_24h_pilot/readout_celltable.csv.gz
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix); library(scuttle)
})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
FIN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_final")
CLN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
OUTD <- file.path(SEN, "processed_data/s1_24h_pilot")
SEED <- 20260924

sm <- read.csv(file.path(SEN, "data/senmayo_saul2022_mouse_final.csv"), stringsAsFactors = FALSE)
SENMAYO_REAL <- unique(trimws(sm$gene))
cat(sprintf("[s1_05a] SenMayo_real = %d genes (saul2022 附表实解)\n", length(SENMAYO_REAL)))

GENES_STAT <- c("Cdkn1a", "Bax", "Mcl1")
FLOOR_MAIN <- as.numeric(read.csv(file.path(FIN, "numi_floors.csv"))$floor_used[1])
FLOORS <- c(lo20 = round(FLOOR_MAIN * 0.8), hi20 = round(FLOOR_MAIN * 1.2))
cat(sprintf("[s1_05a] main floor=%.0f | lo20=%.0f | hi20=%.0f\n", FLOOR_MAIN, FLOORS[1], FLOORS[2]))

metrics <- function(cnt) {                 # cnt: dgCMatrix 基因×细胞 (下采样后)
  lib  <- Matrix::colSums(cnt)
  lnm  <- log1p(t(t(cnt) / pmax(lib, 1)) * 1e4)   # LogNormalize(1e4) 与 x23 逐位同式
  g    <- intersect(SENMAYO_REAL, rownames(cnt))
  sen  <- Matrix::colMeans(lnm[g, , drop = FALSE])
  stat <- do.call(rbind, lapply(GENES_STAT, function(gn)
    if (gn %in% rownames(cnt)) as.numeric(cnt[gn, ]) else rep(NA_real_, ncol(cnt))))
  rownames(stat) <- GENES_STAT
  data.frame(senmayo = as.numeric(sen),
             p21_umis = stat["Cdkn1a", ], bax_umis = stat["Bax", ], mcl1_umis = stat["Mcl1", ],
             nUMI = as.numeric(lib))
}
ds_at <- function(cnt, floor) {            # 确定性下采样 (敏感性版)
  tot  <- Matrix::colSums(cnt)
  prop <- ifelse(tot > floor, floor / tot, 1)
  set.seed(SEED)
  as(scuttle::downsampleMatrix(cnt, prop = prop), "dgCMatrix")
}

ids <- sub("\\.rds$", "", list.files(FIN, pattern = "\\.rds$"))
stopifnot(length(ids) == 28)
cov_n <- 0
out <- data.frame()
for (sid in ids) {
  sceF <- readRDS(file.path(FIN, paste0(sid, ".rds")))
  sceC <- readRDS(file.path(CLN, paste0(sid, ".rds")))
  stopifnot(identical(colnames(sceF), colnames(sceC)))   # 两产物同细胞集
  m_main <- metrics(counts(sceF))
  m_lo   <- metrics(ds_at(counts(sceC), FLOORS[["lo20"]]))
  m_hi   <- metrics(ds_at(counts(sceC), FLOORS[["hi20"]]))
  base <- data.frame(cell_id = paste0(sid, "_", colnames(sceF)),
                     sample = sid, group = sceF$group[1], stringsAsFactors = FALSE)
  out <- rbind(out,
    cbind(base, version = "main", m_main),
    cbind(base, version = "lo20", m_lo),
    cbind(base, version = "hi20", m_hi))
  cov_n <- c(cov_n, length(intersect(SENMAYO_REAL, rownames(sceF)))[1])[-1]
  cat(sprintf("[s1_05a] %-10s cells=%d senmayo基因覆盖=%d\n",
              sid, ncol(sceF), length(intersect(SENMAYO_REAL, rownames(sceF)))))
}
cat(sprintf("[s1_05a] SenMayo 覆盖中位 = %d/%d\n", as.integer(median(cov_n)), length(SENMAYO_REAL)))
write.csv(out, gzfile(file.path(OUTD, "readout_celltable.csv.gz")), row.names = FALSE)
cat(sprintf("[s1_05a] DONE: %d 行 (3版本 × %d 细胞)\n", nrow(out), nrow(out) / 3))
