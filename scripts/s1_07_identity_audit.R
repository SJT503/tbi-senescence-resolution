# s1_07_identity_audit.R — 任务1: 24h 臂细胞身份自查 (marker-based)
# 背景: 2026 bioRxiv 预印本 (AnnoAudit) 报告 GSE269748 星形/OPC 标签在 24h 窗口特异性崩塌
#       Step-1 的 s1_04 标签锚在 GSE269748 fullrun 注释上 → 必须自查 Oligo 身份质量
# 口径: ①marker 阳性率 + 比值; ②双联体/污染份额; ③深度匹配(floor=1921, 多球采样)后复算
# 门禁: 每样本 Oligo 真实细胞数; 不通过不写进手稿
# 用法: Rscript s1_07_identity_audit.R
# 产物: results/s1_identity_audit/{identity_audit_permouse.csv, marker_mean_by_sample.csv, LOG.txt}

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix)
})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
SCED <- file.path(SEN, "processed_data/s1_24h_pilot/sce_clean")
LB   <- file.path(SEN, "processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")
OUTD <- file.path(SEN, "results/s1_identity_audit")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

FLOOR <- 1921L   # 与 s1_03 协议一致
SEED  <- 20260924

L <- character()
log <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = ""); L <<- c(L, s)
}

# ---- 分组 (与 s1_05b 逐字一致) ----
GROUPS <- list(
  naive_M      = c("RNASEQ04","RNASEQ08","RNASEQ10","RNASEQ16","RNASEQ23","RNASEQ24","RNASEQ25"),
  CCI_ipsi_M   = c("RNASEQ03","RNASEQ06","RNASEQ15","RNASEQ17L","RNASEQ18L","RNASEQ19L"),
  CCI_HS_M     = c("RNASEQ01","RNASEQ07","RNASEQ13"),
  rCHI_M       = c("RNASEQ09","RNASEQ11","RNASEQ12"),
  CCI_contra_M = c("RNASEQ17R","RNASEQ18R","RNASEQ19R"),
  CCI_F        = c("RNASEQ20L","RNASEQ21L","RNASEQ22L"),
  naive_F      = c("RNASEQ26","RNASEQ27","RNASEQ28"))
SID2GRP <- setNames(rep(names(GROUPS), lengths(GROUPS)), unlist(GROUPS))
ALL28   <- unlist(GROUPS, use.names = FALSE)

# ---- marker 定义 (Oligo 正 / 髓系负 / 神经元负 / 星形负) ----
M_OLIGO  <- c("Plp1","Mbp","Mag","Mog","Cldn11","Mal")            # 成熟少突
M_MYELO  <- c("C1qa","C1qb","Cx3cr1","P2ry12","Lyz2","Csf1r","Hexb","Ctss")  # 小胶质+巨噬
M_NEUR   <- c("Snap25","Syt1","Rbfox3")
M_ASTRO  <- c("Aqp4","Gfap","Slc1a2")
ANYMARK  <- unique(c(M_OLIGO, M_MYELO, M_NEUR, M_ASTRO, "Cdkn1a"))

# ---- 读入锚标签 ----
lab <- read.csv(gzfile(LB), stringsAsFactors = FALSE)
lab <- lab[, c("cell_id","sample","cell_type")]
log(sprintf("[s1_07] 锚标签 %d 细胞; 类型 %d", nrow(lab), length(unique(lab$cell_type))))

# ---- 多球降采样 (与 s1_03 同法: 按 UMI 比例重分配) ----
# ⚠️ 稀疏累积版: 不可用 matrix(0L, genes, cells) 稠密分配 (2万x2.2万=1.8GB 会 OOM)
downsample_counts <- function(cnt, target, seed) {
  set.seed(seed)
  tot <- Matrix::colSums(cnt)
  keep <- tot > 0
  cnt <- cnt[, keep, drop = FALSE]; tot <- tot[keep]
  newtot <- pmin(tot, target)
  ci <- integer(0); ri <- integer(0); xv <- integer(0)
  for (j in seq_len(ncol(cnt))) {
    idx <- which(cnt[, j] > 0)
    if (!length(idx)) next
    pr <- as.numeric(cnt[idx, j])
    v  <- rmultinom(1L, size = as.integer(round(newtot[j])), prob = pr)[, 1L]
    nz <- v > 0
    if (!any(nz)) next
    ci <- c(ci, rep.int(j, sum(nz))); ri <- c(ri, idx[nz]); xv <- c(xv, v[nz])
  }
  Matrix::sparseMatrix(i = ri, j = ci, x = as.numeric(xv),
                       dims = c(nrow(cnt), ncol(cnt)),
                       dimnames = list(rownames(cnt), colnames(cnt)))
}

summ_block <- function(cnt, cells, tag) {
  # cnt: genes x cells (已与该 cells 对齐); 返回单行汇总
  tot <- Matrix::colSums(cnt); nfeat <- Matrix::colSums(cnt > 0)
  m  <- function(g) if (g %in% rownames(cnt)) as.numeric(cnt[g, ]) else rep(0, ncol(cnt))
  det <- function(g, ...) {
    gs <- c(...); gs <- gs[gs %in% rownames(cnt)]
    if (!length(gs)) return(rep(FALSE, ncol(cnt)))
    Matrix::colSums(cnt[gs, , drop = FALSE] > 0) > 0
  }
  posO <- det(NULL, M_OLIGO); posM <- det(NULL, M_MYELO)
  posN <- det(NULL, M_NEUR);  posA <- det(NULL, M_ASTRO)
  data.frame(
    block = tag, n_cells = ncol(cnt),
    med_umi = median(tot), med_nfeat = median(nfeat),
    # 单 marker 阳性率
    Plp1_pct  = 100 * mean(m("Plp1")  > 0),
    Mbp_pct   = 100 * mean(m("Mbp")   > 0),
    Sox10_pct = if ("Sox10" %in% rownames(cnt)) 100 * mean(m("Sox10") > 0) else NA_real_,
    C1qa_pct  = 100 * mean(m("C1qa")  > 0),
    Cx3cr1_pct= 100 * mean(m("Cx3cr1")> 0),
    Lyz2_pct  = 100 * mean(m("Lyz2")  > 0),
    Snap25_pct= 100 * mean(m("Snap25")> 0),
    Aqp4_pct  = 100 * mean(m("Aqp4")  > 0),
    Cdkn1a_pct= 100 * mean(m("Cdkn1a") > 0),
    # 面板阳性率
    oligo_panel_pct = 100 * mean(posO),
    myeloid_panel_pct = 100 * mean(posM),
    neuron_panel_pct  = 100 * mean(posN),
    astro_panel_pct   = 100 * mean(posA),
    # 身份混杂: 同时阳性
    oligo_and_myeloid_pct = 100 * mean(posO & posM),
    # 均值 (log1p CPM 尺度)
    mean_oligo_panel = mean(colMeans(log1p(t(t(cnt[intersect(M_OLIGO, rownames(cnt)), , drop=FALSE]) /
                                              pmax(tot,1) * 1e4)))),
    mean_myeloid_panel = mean(colMeans(log1p(t(t(cnt[intersect(M_MYELO, rownames(cnt)), , drop=FALSE]) /
                                                pmax(tot,1) * 1e4)))),
    stringsAsFactors = FALSE)
}

permouse <- list(); markmean <- list()

for (sid in ALL28) {
  f <- file.path(SCED, paste0(sid, ".rds"))
  if (!file.exists(f)) { log(sprintf("  !! %s 缺失", sid)); next }
  sce <- readRDS(f)
  cell_id <- paste0(sid, "_", colnames(sce))
  ct <- lab$cell_type[match(cell_id, lab$cell_id)]

  cnt_all <- as(counts(sce), "dgCMatrix")
  genes_ok <- rownames(cnt_all)
  # decontX 污染分数
  dcont <- if ("decontX_contamination" %in% colnames(colData(sce)))
    as.numeric(colData(sce)$decontX_contamination) else rep(NA_real_, ncol(sce))

  # --- 全库(所有已标注细胞) 与 Oligo 子集 ---
  ok_any  <- !is.na(ct) & ct != "Ambiguous"
  sel_oli <- which(ok_any & ct == "Oligo")
  sel_mye <- which(ok_any & ct %in% c("Microglia","Macrophage","LowConf_Microglia"))

  blk <- list()
  if (any(ok_any))  blk[["all_labeled"]] <- cnt_all[, which(ok_any), drop = FALSE]
  if (length(sel_oli)) blk[["Oligo"]]    <- cnt_all[, sel_oli, drop = FALSE]
  if (length(sel_mye)) blk[["Myeloid"]]  <- cnt_all[, sel_mye, drop = FALSE]

  # 深度匹配版 Oligo (floor = 1921)
  if (length(sel_oli) >= 5) {
    co <- cnt_all[, sel_oli, drop = FALSE]
    if (median(Matrix::colSums(co)) > FLOOR) {
      blk[["Oligo_down1921"]] <- downsample_counts(co, FLOOR, SEED + which(ALL28 == sid))
    } else {
      blk[["Oligo_down1921"]] <- co   # 已低于 floor, 不升采样 (与 s1_03 同则)
    }
  }

  for (nm in names(blk)) {
    r <- summ_block(blk[[nm]], NULL, nm)
    r$sample <- sid; r$group <- SID2GRP[[sid]]
    permouse[[length(permouse) + 1]] <- r
  }
  # 污染分数汇总
  if (length(sel_oli)) {
    markmean[[length(markmean) + 1]] <- data.frame(
      sample = sid, group = SID2GRP[[sid]], n_oligo = length(sel_oli),
      med_contam_oligo = median(dcont[sel_oli], na.rm = TRUE),
      p90_contam_oligo = quantile(dcont[sel_oli], .9, na.rm = TRUE),
      med_contam_myeloid = if (length(sel_mye)) median(dcont[sel_mye], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE)
  }
  log(sprintf("  %-10s all=%6d oligo=%5d myeloid=%6d medContam(oligo)=%.4f",
              sid, sum(ok_any), length(sel_oli), length(sel_mye),
              if (length(sel_oli)) median(dcont[sel_oli], na.rm = TRUE) else NA_real_))
  rm(sce, cnt_all); gc(verbose = FALSE)
}

pm <- do.call(rbind, permouse)
mm <- do.call(rbind, markmean)
write.csv(pm, file.path(OUTD, "identity_audit_permouse.csv"), row.names = FALSE)
write.csv(mm, file.path(OUTD, "contamination_by_sample.csv"), row.names = FALSE)

# ---- 汇总: 24h 臂 (CCI_ipsi_M + naive_M) 的 Oligo 身份质量 ----
log("\n================ 汇总: Oligo 身份质量 ================")
oli <- pm[pm$block == "Oligo", ]
for (grp in c("naive_M","CCI_ipsi_M")) {
  s <- oli[oli$group == grp, ]
  if (!nrow(s)) next
  log(sprintf("%-12s n=%d 鼠 | Plp1+ %5.1f%% | Mbp+ %5.1f%% | C1qa+ %5.1f%% | Cx3cr1+ %5.1f%% | Lyz2+ %5.1f%% | Snap25+ %5.1f%%",
              grp, nrow(s), mean(s$Plp1_pct), mean(s$Mbp_pct), mean(s$C1qa_pct),
              mean(s$Cx3cr1_pct), mean(s$Lyz2_pct), mean(s$Snap25_pct)))
  log(sprintf("%-12s   oligo面板+ %5.1f%% | 髓系面板+ %5.1f%% | 双阳 %5.1f%% | 星形面板+ %5.1f%%",
              "", mean(s$oligo_panel_pct), mean(s$myeloid_panel_pct),
              mean(s$oligo_and_myeloid_pct), mean(s$astro_panel_pct)))
}
# 深度匹配前后对比
d1 <- pm[pm$block == "Oligo", ]; d2 <- pm[pm$block == "Oligo_down1921", ]
if (nrow(d2)) {
  log(sprintf("\n深度匹配(1921)前后 (24h 两臂均值):"))
  for (k in c("C1qa_pct","Cx3cr1_pct","Lyz2_pct","myeloid_panel_pct","Cdkn1a_pct")) {
    a <- mean(d1[[k]][d1$group %in% c("naive_M","CCI_ipsi_M")], na.rm = TRUE)
    b <- mean(d2[[k]][d2$group %in% c("naive_M","CCI_ipsi_M")], na.rm = TRUE)
    log(sprintf("  %-20s 匹配前 %6.2f%% → 匹配后 %6.2f%%", k, a, b))
  }
}

writeLines(L, file.path(OUTD, "LOG.txt"))
log("\n[s1_07] DONE → results/s1_identity_audit/")
writeLines(L, file.path(OUTD, "LOG.txt"))
