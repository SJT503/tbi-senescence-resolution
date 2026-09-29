# s1_21c — GSE271769 深度分层判决分析 (攻克 "p21 深度依赖不可判")
# 思路: 检出率 vs nUMI 分箱, 分组画曲线; 组间重叠区内逐箱配对检验 (鼠级 4v4)
#       + 模型验证 (观察低箱 NoTBI vs 解析下采样预测 NoTBI)
# 判决: 曲线贴合 => 差异=深度伪影; D7/D28 系统性抬高 => 真实升高; D1 抬高 => 方法灵敏度自证
suppressWarnings(suppressMessages({}))

SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUTD <- file.path(SEN, "results/s1_gse271769")
LOGF <- file.path(OUTD, "LOG_DEPTHBIN.txt")
log <- function(...) cat(..., "\n", file=LOGF, append=TRUE)
cat("", file=LOGF)

D <- read.csv(gzfile(file.path(OUTD, "celltable.csv.gz")))
stopifnot(nrow(D) == 52478, length(unique(D$sample)) == 16)
log(sprintf("[depthbin] %d 细胞 %d 样本 (门禁通过)", nrow(D), length(unique(D$sample))))

EDGES <- c(0, 600, 1000, 1400, 1800, 2200, 2600, 3000, 3400, Inf)
LABS  <- c("<600","600-1000","1000-1400","1400-1800","1800-2200","2200-2600",
           "2600-3000","3000-3400","3400+")
D$bin <- cut(D$nUMI, breaks=EDGES, labels=LABS, right=FALSE)
MINCELL <- 30   # 每鼠每箱最少细胞数才贡献

# ---- ① 分组×分箱检出率 (pooled + 逐鼠) ----
for (gene in c("Cdkn1a","Cdkn2a")) {
  log(sprintf("\n======== %s: 检出率(%%) vs nUMI 箱, 分组曲线 ========", gene))
  hdr <- sprintf("%-10s", "bin")
  for (g in c("NoTBI","D1","D7","D28")) hdr <- paste0(hdr, sprintf(" %-22s", g))
  log(hdr)
  for (b in LABS) {
    row <- sprintf("%-10s", b)
    for (g in c("NoTBI","D1","D7","D28")) {
      s <- D[D$bin == b & D$time_point == g, ]
      if (nrow(s) >= MINCELL) {
        # 逐鼠均值 -> 组均值 (与主分析同法)
        ag <- tapply(s[[gene]] > 0, s$sample, mean)
        row <- paste0(row, sprintf(" %5.1f%%(n=%5d,%dm)", 100*mean(ag), nrow(s), length(ag)))
      } else row <- paste0(row, sprintf(" %22s", "-"))
    }
    log(row)
  }
}

# ---- ② 重叠区逐箱配对检验 (鼠级, wilcox 4v4) ----
OVERLAP <- c("1000-1400","1400-1800","1800-2200","2200-2600","2600-3000","3000-3400")
for (gene in c("Cdkn1a","Cdkn2a")) {
  log(sprintf("\n======== %s: 重叠区逐箱 D? vs NoTBI (鼠级 wilcox, 每鼠≥%d 细胞) ========", gene, MINCELL))
  for (b in OVERLAP) {
    for (tp in c("D1","D7","D28")) {
      mv <- function(g) {
        s <- D[D$bin == b & D$time_point == g, ]
        ag <- tapply(s[[gene]] > 0, s$sample, function(v) if(length(v)>=MINCELL) mean(v) else NA_real_)
        ag[!is.na(ag)]
      }
      x <- mv(tp); y <- mv("NoTBI")
      if (length(x) >= 3 && length(y) >= 3) {
        p <- tryCatch(wilcox.test(x, y, exact=TRUE)$p.value, error=function(e) NA_real_)
        log(sprintf("  [%s] %s: %.1f%% (n=%dm) vs NoTBI %.1f%% (n=%dm)  差 %+.1fpp  p=%.4g",
                    b, tp, 100*mean(x), length(x), 100*mean(y), length(y),
                    100*(mean(x)-mean(y)), p))
      } else {
        log(sprintf("  [%s] %s: 鼠数不足 (%dv%d), 跳过", b, tp, length(x), length(y)))
      }
    }
  }
}

# ---- ③ 模型验证: 观察 NoTBI 低箱 vs 解析下采样预测 ----
pdet <- function(k, N, F) {
  k <- as.numeric(k); N <- as.numeric(N)
  r <- ifelse(N <= F, as.numeric(k > 0), NA_real_)
  dn <- which(!(N <= F))
  if (length(dn)) {
    kk <- k[dn]; NN <- N[dn]
    p <- ifelse(kk <= 0, 0, NA_real_)
    nz <- kk > 0
    if (any(nz)) {
      lc <- function(n, r2) lgamma(n+1) - lgamma(r2+1) - lgamma(n-r2+1)
      p[nz] <- 1 - exp(lc(NN[nz]-kk[nz], F) - lc(NN[nz], F))
    }
    r[dn] <- p
  }
  pmin(pmax(r, 0), 1)
}
log("\n======== 模型验证: NoTBI 观察低箱 vs 解析下采样 (k/N 交换性假设) ========")
for (F in c(1200, 1600, 2000, 2400, 2800, 3200)) {
  b <- LABS[which(EDGES[-length(EDGES)] <= F & EDGES[-1] > F)[1]]
  obs <- D[D$time_point=="NoTBI" & D$nUMI >= F-200 & D$nUMI < F+200, ]
  ana <- mean(pdet(obs$Cdkn1a, obs$nUMI, F))   # 同一批细胞解析下采样
  log(sprintf("  F=%d: 该窗口 NoTBI 细胞观察检出 %.1f%% (n=%d) | 解析下采样 %.1f%% | 比值 %.2f",
              F, 100*mean(obs$Cdkn1a>0), nrow(obs), 100*ana, mean(obs$Cdkn1a>0)/ana))
}
# 参照: 同窗口 D7/D1/D28 的观察检出 (重叠区曲线的可读版本)
log("\n======== 参照: 各组在 ±200 窗口的观察检出 (与上表 NoTBI 对读) ========")
for (F in c(1200, 1600, 2000, 2400, 2800, 3200)) {
  row <- sprintf("  nUMI~%d:", F)
  for (g in c("NoTBI","D1","D7","D28")) {
    s <- D[D$time_point==g & D$nUMI >= F-200 & D$nUMI < F+200, ]
    row <- paste0(row, if (nrow(s)>=MINCELL)
      sprintf(" %s=%.1f%%(n=%d)", g, 100*mean(s$Cdkn1a>0), nrow(s)) else
      sprintf(" %s=-", g))
  }
  log(row)
}
log("\n[depthbin] DONE")
cat(readLines(LOGF), sep="\n")
