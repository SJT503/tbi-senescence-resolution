# s1_09_anno_audit_bench.R — 用 AnnoAudit 的 D1 判据原样复刻，检查我们的 Oligo 标签
# 目的: 与预印本报的 93.9% self-marker 同尺对比; 并检查成熟 Oligo 是否也有 24h 时点依赖的注释崩坏
# D1 判据(原文): mean log1p(CPM) 打分; 若非目标(髓系)独立基因均值 > 目标(少突)独立基因均值 × 1.5 → 确认污染
#   髓系独立集 = Trem2/Cd68/Lyz2 ; 少突独立集 = Mbp/Plp1/Mog/Mag
# 输出: results/s1_identity_audit/anno_d1_oligo_permouse.csv + LOG
suppressPackageStartupMessages({library(SingleCellExperiment); library(Matrix)})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SCED <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
LB   <- file.path(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")
OUTD <- file.path(SEN, "results/s1_identity_audit")

MYE <- c("Trem2","Cd68","Lyz2")            # 原文 D1 髓系独立集
OLG <- c("Mbp","Plp1","Mog","Mag")         # 原文 D1 少突独立集
TAU <- 1.5

lab <- read.csv(gzfile(LB), stringsAsFactors = FALSE)[, c("cell_id","cell_type")]
sids <- sort(unique(sub("_.*$", "", lab$cell_id)))

L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }
rows <- list()

for (sid in sids) {
  f <- file.path(SCED, paste0(sid, ".rds"))
  if (!file.exists(f)) next
  sce <- readRDS(f)
  cid <- paste0(sid, "_", colnames(sce))
  ct  <- lab$cell_type[match(cid, lab$cell_id)]
  sel <- which(!is.na(ct) & ct == "Oligo")          # 只看被判为 Oligo 的细胞
  if (!length(sel)) { rm(sce); next }
  cnt <- as(counts(sce), "dgCMatrix")[, sel, drop = FALSE]
  tot <- Matrix::colSums(cnt)
  mlcpm <- function(gs) {
    gs <- intersect(gs, rownames(cnt))
    if (!length(gs)) return(rep(NA_real_, ncol(cnt)))
    sub <- as.matrix(cnt[gs, , drop = FALSE])
    cpm <- sweep(sub, 2, pmax(tot, 1), "/") * 1e4
    colMeans(log1p(cpm))
  }
  s_my <- mlcpm(MYE); s_ol <- mlcpm(OLG)
  tau  <- s_my / pmax(s_ol, 1e-9)
  # 自 marker 通过 = 少突自身基因高于髓系 1.5 倍
  self_pass <- tau <= (1 / TAU)
  rows[[sid]] <- data.frame(
    sample = sid, n_oligo = length(sel),
    med_nUMI = median(tot),
    self_oligo_mean = mean(s_ol), myeloid_mean = mean(s_my),
    tau_median = median(tau), tau_p90 = quantile(tau, .9),
    d1_contam_pct = 100 * mean(tau > TAU),          # 确认为髓系污染
    self_pass_pct = 100 * mean(self_pass),          # 自身基因占优
    ctrl = 100 * mean(tau >= (1/TAU) & tau <= TAU), # 中间带(判不了)
    stringsAsFactors = FALSE)
  log(sprintf("  %-10s n=%4d medUMI=%6.0f tau中位=%.3f D1污染=%5.1f%% 自身占优=%5.1f%%",
              sid, length(sel), median(tot), median(tau), 100*mean(tau > TAU), 100*mean(self_pass)))
  rm(sce, cnt); gc(verbose = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUTD, "anno_d1_oligo_permouse.csv"), row.names = FALSE)

log("\n=== AnnoAudit D1 口径汇总 (τ=1.5) ===")
log(sprintf("全部 Oligo (28 样本): 自身占优 %.1f%% | D1确认髓系污染 %.1f%% | 中间带 %.1f%%",
            mean(res$self_pass_pct), mean(res$d1_contam_pct), mean(res$ctrl)))
log("\n对照 AnnoAudit 原文参考值 (其报的成熟 Oligo = 93.9% self-marker, 全时点合并, 未做时点分层):")
log("  OPC 24h: self 16.4% / microglia 70.7% / D1确认 76.2%")
log("  星形 24h: self 14.2% / microglia 56.7% / D1确认 49.2%")
writeLines(L, file.path(OUTD, "LOG_anno_d1.txt"))
cat("\n[s1_09] DONE\n")
