pkgs <- c("Seurat", "celda", "decontX", "scDblFinder", "ggplot2", "Matrix")
for (p in pkgs) {
  cat(p, ":", requireNamespace(p, quietly = TRUE), "\n")
}
if (requireNamespace("Seurat", quietly = TRUE)) cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
