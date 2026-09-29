# s1_10_timecourse_277487.R — GSE277487 时间轴 (Ctrl/6h/2d/4d) 的少突 p21 程序
# 目的: 决定"急性瞬态"故事有没有第二幕 —— 24h 峰之后是否回落
# 数据: results/fullrun/ckpt/ckpt_GSE277487.rds (Seurat, 162,890 细胞, Oligo 12,690)
# 口径(对齐 s1_01→s1_05b, 但本数据集无 raw 液滴 → 省掉 decontX/scDblFinder, 理由见 LOG):
#   ① 复用 ckpt 的既有注释 (本数据集已在 x22 统一重注释)
#   ② Oligo 子集内做去污染敏感性 (decontX 在 12,690 细胞上秒级)
#   ③ nUMI 共同 floor 下采样 → p21 检出率 (鼠级)
#   ④ 统计: 非单调时间轴 → Kruskal-Wallis + 各时点 vs Ctrl 成对 exact MWU (不用 JT)
# 用法: Rscript s1_10_timecourse_277487.R [smoke]
# 产物: results/s1_277487_timecourse/{permouse.csv, LOG.txt}
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(decontX); library(SingleCellExperiment)
})

SEN  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
CKPT <- file.path(SEN, "results/fullrun/ckpt/ckpt_GSE277487.rds")
SMY  <- file.path(SEN, "data/senmayo_saul2022_mouse_final.csv")
OUTD <- file.path(SEN, "results/s1_277487_timecourse"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
smoke <- "smoke" %in% args          # 严格匹配: 只有显式传 smoke 才进 smoke 模式
cat(sprintf("[s1_10] trailingOnly args = [%s] → smoke = %s\n", paste(args, collapse=", "), smoke))
SEED  <- 20260924
L <- character(); log <- function(...) { s <- paste0(...); cat(s,"\n",sep=""); L <<- c(L,s) }

# ---- 读入 ----
obj <- readRDS(CKPT)
log(sprintf("[s1_10] ckpt: %d 细胞 x %d 基因 | 层: %s",
            ncol(obj), nrow(obj), paste(SeuratObject::Layers(obj[["RNA"]]), collapse=",")))
obj$tp <- factor(obj$time_point, levels = c("Ctrl","group1","group2","group3"))
obju <- subset(obj, cells = colnames(obj)[obj$cell_type == "Oligo"])
rm(obj); gc(verbose = FALSE)
log(sprintf("[s1_10] Oligo %d 细胞; 各时点 %s", ncol(obju),
            paste(names(table(obju$tp)), table(obju$tp), sep="=", collapse=" ")))

if (smoke) {
  ks <- unlist(lapply(split(colnames(obju), obju$tp), function(x) head(x, 400)))
  obju <- subset(obju, cells = ks)
  log(sprintf("[s1_10] SMOKE 模式: %d 细胞", ncol(obju)))
}

cnt <- as(SeuratObject::LayerData(obju[["RNA"]], layer = "counts"), "dgCMatrix")

# ---- ② 去污染 (Oligo 子集, 秒级) + 双联体代理检查 ----
MYE <- c("Trem2","Cd68","Lyz2","C1qa","C1qb","Cx3cr1","P2ry12","Csf1r","Hexb")
MYE <- intersect(MYE, rownames(cnt))
log(sprintf("[s1_10] 髓系标记探测 (%d 基因): 检出率 %s", length(MYE),
            paste(sprintf("%s=%.1f%%", MYE, 100*Matrix::rowMeans(cnt[MYE,,drop=FALSE] > 0)), collapse=" ")))

t0 <- Sys.time()
dc <- tryCatch({
  sce <- SingleCellExperiment(list(counts = cnt))
  sce <- decontX(sce, seed = SEED)
  list(contam = as.numeric(sce$decontX_contamination), ok = TRUE)
}, error = function(e) { log(paste("[s1_10] decontX 失败:", conditionMessage(e))); list(contam = rep(NA_real_, ncol(cnt)), ok = FALSE) })
log(sprintf("[s1_10] decontX: ok=%s 中位污染=%.4f (%.1f 分钟)", dc$ok,
            median(dc$contam, na.rm=TRUE), as.numeric(difftime(Sys.time(), t0, units="mins"))))
contam <- dc$contam; names(contam) <- colnames(cnt)

# ---- SenMayo (小鼠符号首字母大写: Cdkn1a; 用 toupper 做大小写无关匹配) ----
sm <- read.csv(SMY, stringsAsFactors = FALSE)
gcol <- grep("gene|symbol|Gene", colnames(sm), ignore.case = TRUE, value = TRUE)[1]
smg_raw <- unique(trimws(as.character(sm[[gcol]])))
smg_raw <- smg_raw[!is.na(smg_raw) & smg_raw != ""]
# 矩阵行名去重（防 duplicate 导致索引错位）
rowU <- toupper(rownames(cnt))
keep_row <- !duplicated(rowU)
if (!all(keep_row)) { log(sprintf("[s1_10] 去重 %d 个重名基因行", sum(!keep_row))); cnt <- cnt[keep_row, ]; rowU <- rowU[keep_row] }
idx <- match(toupper(smg_raw), rowU)
smg <- rownames(cnt)[idx[!is.na(idx)]]
log(sprintf("[s1_10] SenMayo 表 %d 基因 → 矩阵命中 %d", length(smg_raw), length(smg)))
if (length(smg) < 50) log("[s1_10] !! 覆盖过低, 检查基因符号大小写/物种")

# ---- ③ 共同 floor 下采样 ----
tot_raw <- Matrix::colSums(cnt)
med_by_tp <- tapply(tot_raw, obju$tp, median)
FLOOR <- as.integer(floor(min(med_by_tp)))
log(sprintf("[s1_10] 各时点中位 nUMI: %s | 共同 floor = %d",
            paste(names(med_by_tp), as.integer(med_by_tp), sep="=", collapse=" "), FLOOR))

downsample_sparse <- function(m, target, seed) {
  set.seed(seed)
  tot <- Matrix::colSums(m); keep <- tot > 0
  m <- m[, keep, drop=FALSE]; tot <- tot[keep]
  tgt <- pmin(tot, target)
  ci <- integer(0); ri <- integer(0); xv <- integer(0)
  for (j in seq_len(ncol(m))) {
    idx <- which(m[, j] > 0); if (!length(idx)) next
    v <- rmultinom(1L, size = as.integer(round(tgt[j])), prob = as.numeric(m[idx, j]))[, 1L]
    nz <- v > 0; if (!any(nz)) next
    ci <- c(ci, rep.int(j, sum(nz))); ri <- c(ri, idx[nz]); xv <- c(xv, v[nz])
  }
  Matrix::sparseMatrix(i=ri, j=ci, x=as.numeric(xv),
                       dims=c(nrow(m), ncol(m)), dimnames=dimnames(m))
}
cnt_ds <- downsample_sparse(cnt, FLOOR, SEED)
log(sprintf("[s1_10] 下采样后中位 nUMI: %s",
            paste(names(tapply(Matrix::colSums(cnt_ds), obju$tp, median)),
                  as.integer(tapply(Matrix::colSums(cnt_ds), obju$tp, median)), sep="=", collapse=" ")))

# ---- 按版本计算 senmayo / p21 ----
compute <- function(m, tag) {
  tot <- Matrix::colSums(m)
  cpm <- sweep(as.matrix(m[smg, , drop=FALSE]), 2, pmax(tot,1), "/") * 1e4
  sm_score <- unname(colMeans(log1p(cpm)))          # unname: 防 colMeans 带列名 → data.frame 误当行名
  gene_det <- function(g) if (g %in% rownames(m)) as.numeric(m[g, ] > 0) else rep(0, ncol(m))
  cc <- colnames(m)
  grp_v <- as.character(obju$tp)[match(cc, colnames(obju))]
  mouse_v <- as.character(obju$sample)[match(cc, colnames(obju))]
  contam_v <- unname(contam[match(cc, names(contam))])
  data.frame(cell = cc, tp = grp_v, mouse = mouse_v, version = tag,
             nUMI = unname(tot), senmayo = sm_score,
             p21 = gene_det("Cdkn1a"), p16 = gene_det("Cdkn2a"),
             myeloid_any = as.numeric(Matrix::colSums(m[MYE, , drop=FALSE] > 0) > 0),
             contam = contam_v, row.names = NULL, stringsAsFactors = FALSE)
}
df <- do.call(rbind, list(compute(cnt_ds, "down"), compute(cnt, "raw")))
write.csv(df, gzfile(file.path(OUTD, "celltable.csv.gz")), row.names = FALSE)

# ---- ④ 鼠级汇总 + 统计 ----
agg <- function(d, col) {
  # 固定列名为 value: 若用指标名当列名, rbind 时列名不一致 → match.names 报错
  a <- aggregate(d[[col]], by = list(mouse = d$mouse, tp = d$tp), FUN = function(v) {
    if (col %in% c("p21","p16","myeloid_any")) 100*mean(v>0) else mean(v)
  })
  names(a)[3] <- "value"; a
}
res <- NULL
for (v in unique(df$version)) {
  d <- df[df$version == v, ]
  for (cc in c("p21","senmayo","p16","myeloid_any")) {
    a <- agg(d, cc); a$version <- v; a$metric <- cc
    res <- rbind(res, a)
  }
}
res <- res[, c("version","metric","mouse","tp","value")]
write.csv(res, file.path(OUTD, "permouse.csv"), row.names = FALSE)

log(paste0("\n", strrep("=", 72)))
for (v in c("down","raw")) {
  log(sprintf("\n########## version = %s ##########", v))
  for (cc in c("p21","senmayo","p16","myeloid_any")) {
    s <- res[res$version == v & res$metric == cc, ]
    by <- split(s$value, s$tp)
    if (length(unique(s$tp)) < 2) next
    # 每时点均值 + n
    msg <- paste(sprintf("%s: %.3f (n=%d)", names(by), sapply(by, mean), sapply(by, length)), collapse = " | ")
    kw <- tryCatch(kruskal.test(value ~ factor(tp), data = s)$p.value, error = function(e) NA_real_)
    log(sprintf("\n%s:\n   %s", cc, msg))
    log(sprintf("   Kruskal-Wallis p = %.4g", kw))
    for (g in c("group1","group2","group3")) {
      x <- s$value[s$tp == g]; y <- s$value[s$tp == "Ctrl"]
      if (length(x) < 1 || length(y) < 1) next
      p <- tryCatch(wilcox.test(x, y, exact = TRUE)$p.value, error = function(e) NA_real_)
      lab <- c(group1="6h", group2="2d", group3="4d")[g]
      log(sprintf("   vs Ctrl @ %-3s  %6.3f vs %6.3f  p = %.4g  (n=%d v %d)", lab, mean(x), mean(y), p, length(x), length(y)))
    }
  }
}
writeLines(L, file.path(OUTD, "LOG.txt"))
cat("\n[s1_10] DONE\n")
