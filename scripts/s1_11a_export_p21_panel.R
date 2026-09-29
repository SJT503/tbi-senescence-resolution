# s1_11a_export_p21_panel.R — 一次性导出全细胞 p21/p16 小表
# 修正 s1_11 的 I/O 灾难: 原脚本为 11 个细胞型各读一次 6.4GB ckpt (共 11 次)
# 本脚本只读 1 次 ckpt, 导出小表; 下游全部在小表上算
# 输出: results/s1_277487_timecourse/p21_panel_allcells.csv.gz
# 列: cell, cell_type, sample, time_point, nUMI, Cdkn1a, Cdkn2a
suppressPackageStartupMessages({library(Seurat); library(Matrix)})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
CKPT <- file.path(SEN, "results/fullrun/ckpt/ckpt_GSE277487.rds")
OUT  <- file.path(SEN, "results/s1_277487_timecourse/p21_panel_allcells.csv.gz")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
cat(sprintf("[s1_11a] START %s | loading ckpt ...\n", format(t0, "%H:%M:%S")))
obj <- readRDS(CKPT)
cat(sprintf("[s1_11a] loaded: %d cells x %d genes (%.2f min)\n",
            ncol(obj), nrow(obj), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cnt <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
want <- c("Cdkn1a", "Cdkn2a")
have <- intersect(want, rownames(cnt))
cat(sprintf("[s1_11a] 目标基因在场: %s | 缺失: %s\n",
            paste(have, collapse = ", "),
            if (length(setdiff(want, have))) paste(setdiff(want, have), collapse = ", ") else "无"))
stopifnot("Cdkn1a" %in% have)   # p21 是主读出, 缺了必须停

gmat <- as.matrix(cnt[have, , drop = FALSE])

df <- data.frame(
  cell       = colnames(obj),
  cell_type  = as.character(obj$cell_type),
  sample     = as.character(obj$sample),
  time_point = as.character(obj$time_point),
  nUMI       = as.numeric(Matrix::colSums(cnt)),
  stringsAsFactors = FALSE)
for (g in have) df[[g]] <- as.numeric(gmat[g, ])
if (!"Cdkn2a" %in% have) df[["Cdkn2a"]] <- 0L

cat(sprintf("[s1_11a] 导出前自检: 细胞数=%d | Cdkn1a 阳性率=%.2f%% | 中位 nUMI=%.0f\n",
            nrow(df), 100 * mean(df$Cdkn1a > 0), median(df$nUMI)))
stopifnot(nrow(df) == ncol(obj))

write.csv(df, gzfile(OUT), row.names = FALSE)
cat(sprintf("[s1_11a] DONE %s | %d 行 → %s (%.2f min)\n",
            format(Sys.time(), "%H:%M:%S"), nrow(df), OUT,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
