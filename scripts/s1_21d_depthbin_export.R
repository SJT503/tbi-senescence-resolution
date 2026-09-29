# s1_21d — depthbin 结果导出 CSV (供 Biomni 画图; perbin + permouse 双格式)
D <- read.csv(gzfile("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project/results/s1_gse271769/celltable.csv.gz"))
stopifnot(nrow(D) == 52478, length(unique(D$sample)) == 16)

EDGES <- c(0, 600, 1000, 1400, 1800, 2200, 2600, 3000, 3400, Inf)
MIDS  <- c(300, 800, 1200, 1600, 2000, 2400, 2800, 3200, 4000)
D$bin <- cut(D$nUMI, breaks=EDGES, labels=as.character(MIDS), right=FALSE)

rows_bin <- data.frame(); rows_mouse <- data.frame()
for (gene in c("Cdkn1a","Cdkn2a")) {
  D$k <- D[[gene]]
  for (g in c("NoTBI","D1","D7","D28")) {
    for (b in as.character(MIDS)) {
      s <- D[D$bin == b & D$time_point == g, ]
      if (!nrow(s)) next
      ag <- tapply(s$k > 0, s$sample, mean)
      rows_bin <- rbind(rows_bin, data.frame(
        gene=gene, bin_nUMI=as.integer(b), group=g,
        n_cells=nrow(s), n_mice=length(ag),
        det_pct_pooled=round(100*mean(s$k>0), 3),
        det_pct_mousemean=round(100*mean(ag), 3)))
      for (sm in names(ag)) {
        nn <- sum(s$sample == sm)
        rows_mouse <- rbind(rows_mouse, data.frame(
          gene=gene, bin_nUMI=as.integer(b), group=g, sample=sm,
          n_cells=nn, det_pct=round(100*ag[[sm]], 3)))
      }
    }
  }
}
write.csv(rows_bin,   "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project/results/s1_gse271769/depthbin_perbin.csv", row.names=FALSE)
write.csv(rows_mouse, "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project/results/s1_gse271769/depthbin_permouse.csv", row.names=FALSE)
cat(sprintf("perbin %d rows | permouse %d rows\n", nrow(rows_bin), nrow(rows_mouse)))
# 抽查: Cdkn1a bin=1600
chk <- rows_bin[rows_bin$gene=="Cdkn1a" & rows_bin$bin_nUMI==1600, ]
print(chk[, c("group","n_cells","det_pct_pooled","det_pct_mousemean")])
