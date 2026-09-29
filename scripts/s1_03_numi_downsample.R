# s1_03_numi_downsample.R — Step-1 管线第3步: nUMI 组间匹配 (下采样到共同下限)
# 协议§3: "组间中位 nUMI 下采样到共同下限, 敏感性±20%"
# 下限定义(锁死): floor = min(各组 nUMI 中位数) — 用去污染后、去双联体后计数
# 产物: processed_data/s1_24h_pilot/sce_final/<sample>.rds (counts=下采样后) + numi_floors.csv
# 敏感性 ±20% 不在本脚本存矩阵 — s1_05 读出脚本内从 sce_clean 确定性重采样 (seed 固定)
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix); library(scuttle)
})

SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
IN  <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
OUT <- file.path(SEN, "processed_data/s1_24h_pilot/sce_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

smoke <- length(commandArgs(trailingOnly = TRUE)) > 0 && commandArgs(trailingOnly = TRUE)[1] == "smoke"
targets <- if (smoke) c("RNASEQ04","RNASEQ06") else
  sub("\\.rds$", "", list.files(IN, pattern = "\\.rds$"))
stopifnot(smoke || length(targets) == 28)

# ---- 第1遍: 收集每样本每鼠 nUMI 中位数 (smoke 模式用自身中位数演示流程) ----
med_by_sample <- sapply(targets, function(sid) {
  sce <- readRDS(file.path(IN, paste0(sid, ".rds")))
  as.numeric(median(Matrix::colSums(counts(sce))))
})
grp_of <- function(sid) {
  s <- sub("^RNASEQ", "", sid)   # RNASEQ04 -> 04 (正则锚定符撞前缀, 必须先剥前缀)
  if (s %in% c("04","08","10","16","23","24","25")) "naive_M"
  else if (s %in% c("03","06","15","17L","18L","19L")) "CCI_ipsi_M"
  else if (s %in% c("01","07","13")) "CCI_HS_M"
  else if (s %in% c("09","11","12")) "rCHI_M"
  else if (s %in% c("17R","18R","19R")) "CCI_contra_M"
  else if (s %in% c("20L","21L","22L")) "CCI_F"
  else "naive_F"
}
grps <- sapply(names(med_by_sample), grp_of)
grp_med <- tapply(med_by_sample, grps, median)
FLOOR <- if (smoke) round(min(med_by_sample)) else round(min(grp_med))
cat("[s1_03] 组中位 nUMI:\n"); print(round(grp_med))
cat(sprintf("[s1_03] 共同下限 floor = %d UMI/cell%s\n", FLOOR,
            if (smoke) "  (smoke: min of 2 samples, 全量时=min of 7 组中位)" else ""))

# ---- 第2遍: 下采样 + 保存 ----
set.seed(20260924)
for (sid in targets) {
  sce <- readRDS(file.path(IN, paste0(sid, ".rds")))
  cnt <- counts(sce)
  tot <- Matrix::colSums(cnt)
  prop <- ifelse(tot > FLOOR, FLOOR / tot, 1)   # 低于 floor 的细胞保持原样
  ds <- as(scuttle::downsampleMatrix(cnt, prop = prop), "dgCMatrix")  # scuttle 无 downsampleCounts(那是DropletUtils旧名); 保持 double, storage.mode 对 dgCMatrix 报错
  colData(sce)$nUMI_final <- Matrix::colSums(ds)
  assay(sce, "counts") <- as(ds, "dgCMatrix")
  saveRDS(sce, file.path(OUT, paste0(sid, ".rds")))
  cat(sprintf("[s1_03] %-10s medUMI %d -> %d  cells=%d\n", sid,
              as.integer(median(tot)), as.integer(median(Matrix::colSums(ds))), ncol(sce)))
}
write.csv(data.frame(group = names(grp_med), med_nUMI = round(grp_med),
                     floor_used = FLOOR),
          file.path(OUT, if (smoke) "numi_floors_smoke.csv" else "numi_floors.csv"), row.names = FALSE)
cat("[s1_03] DONE:", length(targets), "samples, floor =", FLOOR, "\n")
