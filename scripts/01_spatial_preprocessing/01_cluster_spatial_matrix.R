library(data.table)
library(ggplot2)
library(Matrix)
library(Seurat)

suppressPackageStartupMessages({
})

counts_dir <- "data/processed/Sample00632261/cell_feature_matrix"
meta_csv   <- "data/processed/Sample00632261/cells.csv"
out_dir    <- "data/processed/Sample00632261"

counts   <- readMM(gzfile(file.path(counts_dir, "matrix.mtx.gz")))
features <- fread(file.path(counts_dir, "features.tsv.gz"), sep = "\t", header = FALSE, data.table = FALSE)
barcodes <- fread(file.path(counts_dir, "barcodes.tsv.gz"), sep = "\t", header = FALSE, data.table = FALSE)
counts   <- as(counts, "CsparseMatrix")

gene_names  <- make.unique(as.character(features[[2]]))
barcode_ids <- as.character(barcodes[[1]])
stopifnot(nrow(counts) == length(gene_names), ncol(counts) == length(barcode_ids))
rownames(counts) <- gene_names
colnames(counts) <- barcode_ids

meta <- fread(meta_csv)
stopifnot(all(c("cell_id","x_centroid","y_centroid") %in% names(meta)))
meta$cell_id <- as.character(meta$cell_id)
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_id

keep <- intersect(colnames(counts), rownames(meta))
if (length(keep) < 1000) {
  warning(sprintf("Only %d overlapping cells between counts and metadata. Check IDs.", length(keep)))
}
counts <- counts[, keep, drop = FALSE]
meta   <- meta[keep, , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(meta)))

seu <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA", min.cells = 0, min.features = 0)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")

n_count_per_cell <- Matrix::colSums(counts)
bad_cells <- (!is.finite(n_count_per_cell)) | (n_count_per_cell <= 0)
if (any(bad_cells)) {
  message("Dropping ", sum(bad_cells), " cells with zero/invalid UMI counts.")
  seu <- seu[, !bad_cells]
}
if ("nFeature_RNA" %in% colnames(seu@meta.data)) {
  zero_feat <- seu$nFeature_RNA <= 0 | is.na(seu$nFeature_RNA)
  if (any(zero_feat, na.rm = TRUE)) {
    message("Dropping ", sum(zero_feat, na.rm = TRUE), " cells with zero features.")
    seu <- seu[, !zero_feat]
  }
}
finite_xy <- is.finite(seu$x_centroid) & is.finite(seu$y_centroid)
if (any(!finite_xy)) {
  message("Dropping ", sum(!finite_xy), " cells with non-finite centroids.")
  seu <- seu[, finite_xy]
}

set.seed(1)
use_sct <- TRUE
ok <- TRUE
tryCatch({
  seu <- SCTransform(seu, assay = "RNA", verbose = FALSE)
}, error = function(e) { ok <<- FALSE; message("SCTransform failed: ", conditionMessage(e)) })
if (!ok || !("SCT" %in% Assays(seu))) {
  message("Falling back to NormalizeData pipeline.")
  use_sct <- FALSE
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
}

if (use_sct) {
  seu <- RunPCA(seu, assay = "SCT", verbose = FALSE)
} else {
  seu <- RunPCA(seu, assay = "RNA", verbose = FALSE)
}
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3)
seu <- RunUMAP(seu, dims = 1:30)

p_umap <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + NoLegend()
print(p_umap)


plot_spatial <- function(df, x = "x_centroid", y = "y_centroid", col = NULL,
                         pt_size = 0.35, reverse_y = FALSE, title = NULL) {
  df[[x]] <- as.numeric(df[[x]])
  df[[y]] <- as.numeric(df[[y]])
  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]), , drop = FALSE]

  p <- ggplot(df, aes_string(x = x, y = y, color = col)) +
    geom_point(size = pt_size) +
    coord_fixed() +
    labs(x = "x_centroid", y = "y_centroid", title = title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  if (reverse_y) p <- p + scale_y_reverse()
  p
}

md <- as.data.frame(seu@meta.data)
if ("cell_id" %in% names(md)) names(md)[names(md) == "cell_id"] <- "cell_id_meta"
md$cell_id <- colnames(seu)
md <- md[, !duplicated(names(md)), drop = FALSE]

p_spatial_clusters <- plot_spatial(md, col = "00632261 (Control) seurat_clusters", pt_size = 0.3, reverse_y = FALSE, title = "00632261 (Control) Cluster")
print(p_spatial_clusters)

feature_gene <- "MALAT1"
if (feature_gene %in% rownames(seu)) {
  expr <- FetchData(seu, vars = feature_gene, slot = "data")
  md$feature_expr <- expr[[1]]
  p_spatial_feature <- plot_spatial(md, col = "feature_expr", pt_size = 0.3, reverse_y = FALSE, title = feature_gene) +
    labs(color = feature_gene)
  print(p_spatial_feature)
} else {
  message(sprintf("Gene '%s' not found; skipping feature map.", feature_gene))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(seu, file.path(out_dir, "seurat_spatial_like.rds"))
ggsave(file.path(out_dir, "umap_clusters.png"), plot = p_umap, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "spatial_clusters.png"), plot = p_spatial_clusters, width = 6, height = 6, dpi = 300)
