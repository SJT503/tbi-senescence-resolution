# s1_26 — GSE271769 p16+/p16- 表型复验 (Fig3E 的第二数据集复验; 复刻 s1_21 清洗 + 扩展基因)
suppressPackageStartupMessages({library(SingleCellExperiment); library(decontX); library(scDblFinder)})
SEN <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/senescence_project"
OUT <- file.path(SEN, "results/s1_atlas_gapfill")
LOG <- file.path(OUT, "LOG_271769_PHENO.txt"); log <- function(...) cat(...,"\n",file=LOG,append=TRUE)
cat("",file=LOG); SEED <- 20260925

CKPT <- file.path(SEN, "results/fullrun/ckpt/ckpt_GSE271769.rds")
obj <- readRDS(CKPT)
cnt_all <- SeuratObject::LayerData(obj[["RNA"]], layer="counts")
mg_mask <- obj$cell_type == "Microglia"
log(sprintf("[s1_26] %d 细胞 | MG: %d", ncol(obj), sum(mg_mask)))
hexb <- 100*mean(as.numeric(cnt_all["Hexb", mg_mask] > 0)); snap25 <- 100*mean(as.numeric(cnt_all["Snap25", mg_mask] > 0))
log(sprintf("[s1_26] 身份自检: Hexb+ %.1f%% Snap25+ %.1f%%", hexb, snap25))
stopifnot(hexb > 80, snap25 < 10)

PANEL <- c("Cdkn1a","Cdkn2a","Cst7","Itgax","Cd9","Trem2","Apoe","Lpl","Spp1","P2ry12",
           "Tmem119","Stat1","Irf7","Bcl2","Mcl1","Bax","Lyz2","C1qa","Hexb")
smps <- unique(as.character(obj$sample[mg_mask]))
parts <- list()
for (sm in smps) {
  cells_in <- which(mg_mask & as.character(obj$sample) == sm)
  m <- cnt_all[intersect(rownames(cnt_all), PANEL), cells_in, drop=FALSE]
  nz <- Matrix::colSums(cnt_all[, cells_in]) > 0
  m <- m[, nz, drop=FALSE]
  if (ncol(m) < 50) next
  tp <- as.character(obj$time_point[cells_in])[nz][1]
  sce <- SingleCellExperiment(list(counts = cnt_all[, cells_in][, nz]))
  sce <- decontX(sce, seed = SEED)
  cln <- round(as(assay(sce, "decontXcounts"), "dgCMatrix"))
  cln@x[cln@x < 0] <- 0
  sce2 <- scDblFinder(SingleCellExperiment(list(counts = cln)))
  singlet <- colData(sce2)$scDblFinder.class == "singlet"
  sub <- as.matrix(cln[intersect(rownames(cln), PANEL), singlet, drop=FALSE])
  df <- data.frame(cell_id = colnames(sub), sample = sm, time_point = tp,
                   nUMI = as.numeric(cln[, singlet] |> Matrix::colSums()), t(sub > 0))
  parts[[sm]] <- df
  log(sprintf("  %-22s n=%5d singlet | Cdkn2a+ %.2f%%", sm, nrow(df), 100*mean(df$Cdkn2a)))
}
D <- do.call(rbind, parts)
write.csv(D, gzfile(file.path(OUT, "gse271769_mg_phenotype_cells.csv.gz")), row.names=FALSE)

# p16+ vs p16- 检出率 (D7/D28; 描述性 + 深度比)
for (tp in c("D7","D28")) {
  s <- D[D$time_point == tp, ]
  pos <- s[s$Cdkn2a, ]; neg <- s[!s$Cdkn2a, ]
  log(sprintf("\n[%s] p16+ n=%d vs p16- n=%d | 中位 nUMI %.0f vs %.0f (比值 %.2f)",
              tp, nrow(pos), nrow(neg), median(pos$nUMI), median(neg$nUMI), median(pos$nUMI)/median(neg$nUMI)))
  for (g in setdiff(PANEL, c("Cdkn2a","Hexb"))) {
    log(sprintf("    %-8s %6.1f%% vs %6.1f%%  (%+.1f pp)", g, 100*mean(pos[[g]]), 100*mean(neg[[g]]),
                100*(mean(pos[[g]])-mean(neg[[g]]))))
  }
}
log("\n[s1_26] DONE")
cat(readLines(LOG), sep="\n")
