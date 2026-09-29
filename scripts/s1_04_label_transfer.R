# s1_04_label_transfer.R — Step-1 管线第4步: 细胞类型标注 (协议§3)
#   13 个已用样本: fullrun 锚标签直连 (join by cell_id = sample_barcode)
#   15 个未用样本: Seurat 标签转移 (reference=13样本锚, k.weight=30, 置信度<0.5→Ambiguous)
# 门禁(协议§3): 转移后各型比例与 fullrun 版差异 >2 倍 → stop 排查 (仅对13直连样本可查)
# 用法: Rscript s1_04_label_transfer.R [smoke]   smoke = 转移只用 09(rCHI)+26(naive_F) 2样本
# 产物: processed_data/s1_24h_pilot/labels/{all_cells_labels.csv.gz, prop_gate_check.csv, confidence_by_group.csv}

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(SingleCellExperiment)
})

SEN   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
IN    <- file.path(SEN, "processed_data/s1_24h_pilot/sce_final")
OUT   <- file.path(SEN, "processed_data/s1_24h_pilot/labels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED  <- 20260924

# ---- 样本角色表 (锚覆盖实查 2026-09-24: 13 直连 = naive_M 7 + CCI_ipsi_M 6) ----
DIRECT <- c("RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25",
            "RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L")
TRANSFER <- c("RNASEQ01","RNASEQ07","RNASEQ13","RNASEQ09","RNASEQ11","RNASEQ12",
              "RNASEQ17R","RNASEQ18R","RNASEQ19R","RNASEQ20L","RNASEQ21L","RNASEQ22L",
              "RNASEQ26","RNASEQ27","RNASEQ28")
smoke <- length(commandArgs(trailingOnly = TRUE)) > 0 && commandArgs(trailingOnly = TRUE)[1] == "smoke"
query_ids <- if (smoke) c("RNASEQ09","RNASEQ26") else TRANSFER
all_ids <- c(DIRECT, TRANSFER)
stopifnot(smoke || length(list.files(IN, pattern = "\\.rds$")) == 28)

# ---- 读入 sce_final, 统一 cell_id = sample_barcode (与锚格式一致) ----
load_sce <- function(sid) {
  sce <- readRDS(file.path(IN, paste0(sid, ".rds")))
  colnames(sce) <- paste0(sid, "_", colnames(sce))
  sce
}
cat("[s1_04] loading 13 direct +", length(query_ids), "query samples ...\n")
sces <- c(lapply(DIRECT, load_sce), lapply(query_ids, load_sce))
names(sces) <- c(DIRECT, query_ids)

# ---- 锚标签 (fullrun, 只取 13 直连样本的细胞) ----
anch <- read.csv(gzfile(file.path(SEN, "results/fullrun/x23_scores_GSE269748.csv.gz")),
                 stringsAsFactors = FALSE,
                 colClasses = c(cell_id = "character", cell_type = "character"))
anch <- anch[, c("cell_id", "cell_type")]
anch <- anch[grepl("^RNASEQ(04|08|10|16|23|24|25|03|06|15|17L|18L|19L)_", anch$cell_id), ]
stopifnot(nrow(anch) > 0)

# ---- 直接标注 13 样本 + 比例门禁 ----
mk_df <- function(sid, type_vec, method, conf) {
  data.frame(cell_id = colnames(sces[[sid]]), sample = sid,
             group = sces[[sid]]$group[1],
             cell_type = type_vec, method = method, confidence = conf,
             stringsAsFactors = FALSE)
}
direct_parts <- lapply(DIRECT, function(sid) {
  hit <- match(colnames(sces[[sid]]), anch$cell_id)
  typ <- ifelse(is.na(hit), NA_character_, anch$cell_type[hit])
  data.frame(cell_id = colnames(sces[[sid]]), sample = sid,
             group = sces[[sid]]$group[1], cell_type = typ,
             method = "direct", confidence = 1.0, stringsAsFactors = FALSE)
})
direct_df <- do.call(rbind, direct_parts)
miss_rate <- mean(is.na(direct_df$cell_type))
cat(sprintf("[s1_04] direct join: %d cells, 锚未覆盖(QC/双联体已删) = %.1f%%\n",
            nrow(direct_df), 100 * miss_rate))
if (miss_rate > 0.5) stop("直连缺失率>50%: cell_id 格式或锚版本异常, 排查后再跑")

# 比例门禁: 同样本 fullrun 锚比例 vs 管线后直连比例, 任一 ≥1% 型 max/min > 2 → stop
gate_rows <- do.call(rbind, lapply(DIRECT, function(sid) {
  a <- table(anch$cell_type[grepl(paste0("^", sid, "_"), anch$cell_id)]); a <- a / sum(a)
  d <- table(direct_df$cell_type[direct_df$sample == sid], useNA = "no"); d <- d / sum(d)
  tp <- union(names(a), names(d))
  data.frame(sample = sid, cell_type = tp,
             prop_fullrun = as.numeric(a[tp]), prop_pipeline = as.numeric(d[tp]))
}))
gate_rows$ratio <- pmax(gate_rows$prop_fullrun, gate_rows$prop_pipeline) /
  pmax(gate_rows$prop_fullrun, gate_rows$prop_pipeline, 0.005)
viol <- gate_rows[gate_rows$prop_fullrun >= 0.01 | gate_rows$prop_pipeline >= 0.01, ]
viol <- viol[viol$ratio > 2, ]
cat(sprintf("[s1_04] 比例门禁: %d 项 >2倍\n", nrow(viol)))
if (nrow(viol) > 0) { print(viol); stop("比例门禁 FAIL: 转移后比例与 fullrun 差>2倍, 停止排查") }

# ---- Seurat 标签转移 (reference = 13 直连样本当前管线细胞 + 锚标签) ----
set.seed(SEED)
ref_cells <- direct_df[!is.na(direct_df$cell_type), ]
mat_list <- lapply(DIRECT, function(sid) counts(sces[[sid]]))
genes <- rownames(mat_list[[1]])
mat_ref <- do.call(cbind, mat_list)
stopifnot(all(sapply(mat_list, function(m) identical(rownames(m), genes))))
ref <- CreateSeuratObject(counts = mat_ref, meta.data =
        data.frame(cell_type = ref_cells$cell_type[match(colnames(mat_ref), ref_cells$cell_id)],
                   row.names = colnames(mat_ref)))
ref <- subset(ref, cells = rownames(ref@meta.data)[!is.na(ref@meta.data$cell_type)])
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)  # 默认只 scale 2000 HVG — 全基因会稠密化 ~15GB OOM
ref <- RunPCA(ref, npcs = 30, verbose = FALSE)

mat_qry <- do.call(cbind, lapply(query_ids, function(sid) counts(sces[[sid]])))
qry <- CreateSeuratObject(counts = mat_qry)
qry <- NormalizeData(qry, verbose = FALSE)
cat(sprintf("[s1_04] transfer: ref=%d cells (%d types), query=%d cells (%d samples)\n",
            ncol(ref), length(unique(ref$cell_type)), ncol(qry), length(query_ids)))
tc0 <- Sys.time()
anchors <- FindTransferAnchors(reference = ref, query = qry, dims = 1:30,
                               reference.reduction = "pca", verbose = TRUE)
pred <- TransferData(anchors, refdata = ref$cell_type, dims = 1:30,
                     k.weight = 30, verbose = FALSE)  # 协议: k=30; weight.reduction 用默认 pcaproject(锚内已投影)
cat(sprintf("[s1_04] transfer done in %.1f min\n", as.numeric(difftime(Sys.time(), tc0, units = "mins"))))

conf <- pred$prediction.score.max
typ  <- ifelse(conf < 0.5, "Ambiguous", pred$predicted.id)  # 协议: 置信度<0.5 → Ambiguous
stopifnot(length(typ) == ncol(qry), length(conf) == ncol(qry))
qry_df <- data.frame(cell_id = colnames(qry),
                     sample = sub("_.*$", "", colnames(qry)),
                     group = rep(vapply(query_ids, function(s) sces[[s]]$group[1], ""),
                                 times = vapply(query_ids, function(s) ncol(sces[[s]]), 1L)),
                     cell_type = typ, method = "transfer", confidence = round(conf, 4),
                     stringsAsFactors = FALSE)
stopifnot(all(qry_df$sample == sub("_.*$", "", qry_df$cell_id)))   # group 对位防线
cat("[s1_04] 转移置信度分布(按组):\n")
print(aggregate(confidence ~ group, data = qry_df, FUN = function(x)
  round(c(med = median(x), p25 = quantile(x, .25), ambig_pct = 100 * mean(x < 0.5)), 3)))

# ---- 汇总输出 ----
all_labels <- rbind(
  direct_df[!is.na(direct_df$cell_type), ],
  qry_df)
write.csv(all_labels, gzfile(file.path(OUT, if (smoke) "all_cells_labels_smoke.csv.gz"
                                              else "all_cells_labels.csv.gz")), row.names = FALSE)
write.csv(gate_rows, file.path(OUT, if (smoke) "prop_gate_check_smoke.csv"
                                          else "prop_gate_check.csv"), row.names = FALSE)
write.csv(aggregate(confidence ~ sample + group, data = qry_df, FUN = function(x)
  round(c(med_conf = median(x), pct_lt05 = 100 * mean(x < 0.5)), 3)),
  file.path(OUT, if (smoke) "confidence_by_group_smoke.csv" else "confidence_by_group.csv"),
  row.names = FALSE)
cat(sprintf("[s1_04] DONE: %d labeled cells (direct %d + transfer %d)\n",
            nrow(all_labels), sum(all_labels$method == "direct"), sum(all_labels$method == "transfer")))
print(table(all_labels$cell_type))
