# s1_23 — MG SASP + 细胞周期评分提取 (24h 批, p21+/p21- 对比; 供新补充图)
suppressPackageStartupMessages({library(SingleCellExperiment)})
SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUT <- file.path(SEN, "results/s1_atlas_gapfill"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
LOG <- file.path(OUT, "LOG_SASP.txt"); log <- function(...) cat(...,"\n",file=LOG,append=TRUE)
cat("",file=LOG)

lab <- read.csv(gzfile(file.path(SEN,"processed_data/s1_24h_pilot/labels/all_cells_labels.csv.gz")))
mg <- lab[lab$cell_type=="Microglia", c("cell_id","sample","group")]
log(sprintf("[s1_23] MG 标签: %d 细胞 %d 样本", nrow(mg), length(unique(mg$sample))))

SASP <- c("Il6","Il1a","Il1b","Tnf","Ccl2","Ccl3","Ccl4","Ccl7","Ccl12","Cxcl1","Cxcl2","Cxcl10",
          "Mmp3","Mmp9","Mmp13","Serpine1","Plau","Igfbp3","Igfbp5","Igfbp7","Ptgs2","Vegfa",
          "Gdf15","Tgfb1","Cxcl12","Csf1","Mif")
SG   <- c("Mcm2","Mcm4","Mcm5","Mcm6","Pcna","Tyms","Rrm1","Rrm2","Ung","Ube2t","Fen1","Prim1")
G2M  <- c("Ccnb1","Ccnb2","Cdk1","Top2a","Bub1","Bub1b","Aurkb","Cenpf","Cenpa","Nusap1","Ube2c","Plk1")

parts <- list()
for (sm in sort(unique(mg$sample))) {
  f <- file.path(SEN, sprintf("processed_data/s1_24h_pilot/sce_clean/%s.rds", sm))
  sce <- readRDS(f)
  cm <- as.matrix(counts(sce)[intersect(rownames(sce), c(SASP,SG,G2M,"Cdkn1a","Cdkn2a","Hexb")), , drop=FALSE])
  keep <- paste0(sm, "_", colnames(sce)) %in% mg$cell_id   # sce 条码无样本前缀, labels 有
  cm <- cm[, keep, drop=FALSE]
  colnames(cm) <- paste0(sm, "_", colnames(cm))
  parts[[sm]] <- cbind(t(cm), nUMI_total = sce$nUMI_raw[keep])
  rm(sce); gc(verbose=FALSE)
}
X <- do.call(rbind, parts)
idx <- match(rownames(X), mg$cell_id)
df <- data.frame(cell_id=rownames(X), sample=mg$sample[idx], group=mg$group[idx],
                 as.data.frame(X), check.names=FALSE)
df$p21 <- df$Cdkn1a > 0; df$p16 <- df$Cdkn2a > 0
genes_found <- intersect(c(SASP,SG,G2M), colnames(df))
log(sprintf("[s1_23] 基因命中: SASP %d/%d | S %d/%d | G2M %d/%d",
             sum(intersect(SASP,colnames(df)) %in% colnames(df)), length(SASP),
             length(intersect(SG,colnames(df))), length(SG), length(intersect(G2M,colnames(df))), length(G2M)))
# 模块评分 = 每细胞基因平均检出 (binary mean)
mod <- function(df, genes) apply(as.matrix(df[, intersect(genes, colnames(df))]>0), 1, mean)
df$SASP_score <- mod(df, SASP); df$S_score <- mod(df, SG); df$G2M_score <- mod(df, G2M)
keepcols <- intersect(c("cell_id","sample","group","p21","p16","SASP_score","S_score","G2M_score","nUMI_total",SASP,SG,G2M), colnames(df))
log(sprintf("[s1_23] 保留 %d 列; 缺失: %s", length(keepcols),
            paste(setdiff(c("cell_id","sample","group","p21","p16","SASP_score","S_score","G2M_score"), colnames(df)), collapse=",")))
write.csv(df[, keepcols], gzfile(file.path(OUT,"mg24_sasp_cellcycle.csv.gz")), row.names=FALSE)

# 逐鼠汇总: p21+ vs p21- 的 SASP/cellcycle
agg <- aggregate(cbind(SASP_score,S_score,G2M_score) ~ sample + p21, data=df, FUN=mean)
write.csv(agg, file.path(OUT,"mg24_sasp_cc_permouse.csv"), row.names=FALSE)
for (v in c("SASP_score","S_score","G2M_score")) {
  w <- reshape(agg, idvar="sample", timevar="p21", direction="wide")
  pos <- w[[paste0(v,".TRUE")]]; neg <- w[[paste0(v,".FALSE")]]
  log(sprintf("  [p21±] %s: p21+ %.4f vs p21- %.4f | 差 %+.4f | 方向一致 %d/%d",
              v, mean(pos), mean(neg), mean(pos-neg), sum(pos>neg), length(pos)))
}
log("[s1_23] DONE")
cat(readLines(LOG), sep="\n")
