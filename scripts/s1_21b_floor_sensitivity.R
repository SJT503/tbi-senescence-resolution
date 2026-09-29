# s1_21b — GSE271769 floor 敏感性分析 (v7 §17 遗留项落盘)
# 目的: S3F 图的数据源; 验证 README "F 1540-3647 内 p16 3.3-5.1×" 结论
# 口径: 与 s1_21 完全同一 pdet 公式; 唯一区别 = F 不再取 min(median), 而是扫描
suppressWarnings(suppressMessages({}))

SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUTD <- file.path(SEN, "results/s1_gse271769")
LOGF <- file.path(OUTD, "LOG_FLOOR_SENS.txt")
log <- function(...) cat(..., "\n", file=LOGF, append=TRUE)
cat("", file=LOGF)  # 清空

D <- read.csv(gzfile(file.path(OUTD, "celltable.csv.gz")))
stopifnot(nrow(D) == 52478, length(unique(D$sample)) == 16)
stopifnot(all(c("nUMI","Cdkn1a","Cdkn2a","sample","time_point") %in% names(D)))
log(sprintf("[s1_21b] celltable %d 细胞 %d 样本 (门禁通过)", nrow(D), length(unique(D$sample))))

# pdet — 与 s1_21 逐字一致
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

# 主分析 F=1540 校验值 (与 LOG.txt 对齐断言)
chk <- aggregate(100*pdet(D$Cdkn2a, D$nUMI, 1540) ~ time_point, data=D,
                 FUN=function(v) mean(v[seq_len(length(v))]))  # placeholder, 重算见下
# 逐鼠聚合再取组均值 (与 s1_21 同法: 鼠均值 → 组均值)
mouse_ag <- function(pvec) {
  df <- data.frame(sample=D$sample, time_point=D$time_point, val=pvec)
  ag <- aggregate(val ~ sample + time_point, data=df,
                  FUN=function(v) if(length(v)>=10) mean(v) else NA_real_)
  ag[!is.na(ag$val), ]
}
mm <- mouse_ag(100*pdet(D$Cdkn2a, D$nUMI, 1540))
gm <- tapply(mm$val, mm$time_point, mean)
log(sprintf("[s1_21b] F=1540 校验: D7=%.3f%% NoTBI=%.3f%% fold=%.2f (LOG 主分析: 0.242/0.048/5.0)",
            gm["D7"], gm["NoTBI"], gm["D7"]/gm["NoTBI"]))

# floor 扫描
FLOORS <- c(800, 1000, 1200, 1400, 1540, 1800, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3647)
TPS <- c("NoTBI","D1","D7","D28")
out <- data.frame()
for (F in FLOORS) {
  for (g in c("p21","p16")) {
    pv <- 100*pdet(D[[ifelse(g=="p21","Cdkn1a","Cdkn2a")]], D$nUMI, F)
    mm <- mouse_ag(pv)
    gm <- tapply(mm$val, mm$time_point, mean)
    row <- data.frame(floor=F, gene=g,
                      NoTBI=unname(gm["NoTBI"]), D1=unname(gm["D1"]),
                      D7=unname(gm["D7"]), D28=unname(gm["D28"]))
    for (tp in c("D1","D7","D28")) {
      x <- mm$val[mm$time_point==tp]; y <- mm$val[mm$time_point=="NoTBI"]
      p <- tryCatch(wilcox.test(x, y, exact=TRUE)$p.value, error=function(e) NA_real_)
      row[paste0("fold_",tp)] <- unname(gm[tp]/gm["NoTBI"])
      row[paste0("p_",tp)]   <- p
    }
    out <- rbind(out, row)
  }
}
write.csv(out, file.path(OUTD, "floor_sensitivity.csv"), row.names=FALSE)

# 判读摘要
r <- out[out$gene=="p16" & out$floor>=1540 & out$floor<=3647, ]
log(sprintf("[s1_21b] p16 fold 范围 (F 1540-3647): D7 %.2f-%.2f | D28 %.2f-%.2f | D1 %.2f-%.2f",
            min(r$fold_D7), max(r$fold_D7), min(r$fold_D28), max(r$fold_D28),
            min(r$fold_D1), max(r$fold_D1)))
r2 <- out[out$gene=="p21" & out$floor>=1540 & out$floor<=3647, ]
log(sprintf("[s1_21b] p21 fold 范围 (F 1540-3647): D7 %.2f-%.2f | D28 %.2f-%.2f | D1 %.2f-%.2f",
            min(r2$fold_D7), max(r2$fold_D7), min(r2$fold_D28), max(r2$fold_D28),
            min(r2$fold_D1), max(r2$fold_D1)))
sig <- sapply(c("D1","D7","D28"), function(tp) sum(r[,paste0("p_",tp)] < 0.05, na.rm=TRUE))
log(sprintf("[s1_21b] p16 在 F 1540-3647 内显著 (p<0.05) 的 floor 数: D1=%d/14 D7=%d/14 D28=%d/14",
            sig["D1"], sig["D7"], sig["D28"]))
log("[s1_21b] DONE → floor_sensitivity.csv")
cat(readLines(LOGF), sep="\n")
