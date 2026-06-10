library(Matrix); library(data.table); library(SpatialExperiment); library(scater); library(scran); library(Seurat); library(ggplot2)

counts_dir <- "data/processed/cell_feature_matrix"
meta_csv   <- "data/processed/cells.csv"

counts <- readMM(gzfile(file.path(counts_dir, "matrix.mtx.gz")))
features <- fread(file.path(counts_dir, "features.tsv.gz"), sep = "\t", header = FALSE, data.table = FALSE)
barcodes <- fread(file.path(counts_dir, "barcodes.tsv.gz"), sep = "\t", header = FALSE, data.table = FALSE)

gene_names <- make.unique(as.character(features[[2]]))  # V2 are gene symbols
barcode_ids <- as.character(barcodes[[1]])
rownames(counts) <- gene_names
colnames(counts) <- barcode_ids

stopifnot(nrow(counts) == nrow(features))
stopifnot(ncol(counts) == nrow(barcodes))

meta <- fread(meta_csv) |> as.data.frame()
stopifnot("cell_id" %in% names(meta), "x_centroid" %in% names(meta), "y_centroid" %in% names(meta))
meta$cell_id <- as.character(meta$cell_id)
rownames(meta) <- meta$cell_id

keep <- intersect(colnames(counts), rownames(meta))
if (length(keep) < 1000) {
  warning(sprintf("Only %d overlapping cells between counts and metadata. Double-check IDs.", length(keep)))
}
counts <- counts[, keep, drop = FALSE]
meta   <- meta[keep, , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(meta)))

spe <- SpatialExperiment(
  assays  = list(counts = counts),
  colData = meta
)
spatialCoords(spe) <- as.matrix(meta[, c("x_centroid", "y_centroid")])

set.seed(1)
spe <- logNormCounts(spe)
spe <- runPCA(spe)
g <- scran::buildSNNGraph(spe, use.dimred = "PCA", k = 15)
cl <- igraph::cluster_walktrap(g)$membership
colLabels(spe) <- factor(cl)

df_plot <- cbind(as.data.frame(spatialCoords(spe)), cluster = colLabels(spe))
p_spatial <- ggplot(df_plot, aes(x = x_centroid, y = y_centroid, color = cluster)) +
  geom_point(size = 0.4) + coord_fixed() + theme_minimal()
print(p_spatial)

seu <- CreateSeuratObject(counts = counts, meta.data = meta)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")
seu <- NormalizeData(seu) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
seu <- FindNeighbors(seu, dims = 1:30) |> FindClusters(resolution = 0.6) |> RunUMAP(dims = 1:30)
DimPlot(seu, group.by = "seurat_clusters", label = TRUE) + NoLegend()

saveRDS(spe, "data/processed/spe_aligned.rds")
saveRDS(seu, "data/processed/seurat_aligned.rds")
