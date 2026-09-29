# s1_21c2 — depthbin 补充: ①深细胞单独下采样 ②3400+箱正式检验
D <- read.csv(gzfile("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project/results/s1_gse271769/celltable.csv.gz"))
pdet <- function(k, N, F) {
  r <- ifelse(N <= F, as.numeric(k > 0), NA_real_)
  dn <- which(!(N <= F))
  if (length(dn)) {
    kk <- k[dn]; NN <- N[dn]; p <- rep(0, length(dn)); nz <- kk > 0
    lc <- function(n, r2) lgamma(n+1) - lgamma(r2+1) - lgamma(n-r2+1)
    p[nz] <- 1 - exp(lc(NN[nz]-kk[nz], F) - lc(NN[nz], F)); r[dn] <- p
  }
  pmin(pmax(r, 0), 1)
}
L <- file("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project/results/s1_gse271769/LOG_DEPTHBIN.txt", open="a")
cat("\n======== 补1: 各组深(N>=3400)单独下采样到F=1540 vs 该组浅箱(1400-1800)观察 ========\n", file=L)
for (g in c("NoTBI","D1","D7","D28")) {
  deep <- D[D$time_point==g & D$nUMI>=3400, ]
  shal <- D[D$time_point==g & D$nUMI>=1400 & D$nUMI<1800, ]
  cat(sprintf("  %-6s 深(n=%4d)下采样到1540 = %5.1f%% | 浅箱观察(n=%4d) = %5.1f%%\n",
              g, nrow(deep), 100*mean(pdet(deep$Cdkn1a, deep$nUMI, 1540)),
              nrow(shal), 100*mean(shal$Cdkn1a>0)), file=L)
}
cat("\n======== 补2: 3400+ 箱 D? vs NoTBI (鼠级 wilcox, 每鼠>=30 细胞) ========\n", file=L)
s <- D[D$nUMI>=3400, ]
for (tp in c("D1","D7","D28")) {
  mv <- function(g) { x <- s[s$time_point==g, ]
    ag <- tapply(x$Cdkn1a>0, x$sample, function(v) if(length(v)>=30) mean(v) else NA_real_)
    ag[!is.na(ag)] }
  x <- mv(tp); y <- mv("NoTBI")
  p <- tryCatch(wilcox.test(x, y, exact=TRUE)$p.value, error=function(e) NA_real_)
  cat(sprintf("  [3400+] %s: %.1f%% (n=%dm) vs NoTBI %.1f%% (n=%dm)  差 %+.1fpp  p=%.4g\n",
              tp, 100*mean(x), length(x), 100*mean(y), length(y),
              100*(mean(x)-mean(y)), p), file=L)
}
close(L)
cat("done\n")
