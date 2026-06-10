library(dplyr)
library(Matrix)
library(reticulate)
library(Seurat)

setwd("data/processed/Sample00632301")

counts <- readMM("cell_feature_matrix/matrix.mtx.gz")

features <- read.delim("cell_feature_matrix/features.tsv.gz",
                       header = FALSE, stringsAsFactors = FALSE)
colnames(features) <- c("feature_id", "gene_name", "feature_type")

barcodes <- read.delim("cell_feature_matrix/barcodes.tsv.gz",
                       header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

rownames(counts) <- features$feature_id
colnames(counts) <- barcodes$barcode


zarr <- import("zarr")
cells_store <- zarr$open("cells.zarr")

cell_ids <- py_to_r(cells_store$attrs$get("cell_id", NULL))
if (is.null(cell_ids)) {
  cell_ids <- py_to_r(cells_store[["obs"]][["cell_id"]])
}

x_coord <- py_to_r(cells_store[["obs"]][["x_centroid"]])
y_coord <- py_to_r(cells_store[["obs"]][["y_centroid"]])

cells <- data.frame(
  barcode = cell_ids,
  x_centroid = x_coord,
  y_centroid = y_coord,
  stringsAsFactors = FALSE
)

cells <- cells %>%
  filter(barcode %in% colnames(counts)) %>%
  distinct(barcode, .keep_all = TRUE)

cells <- cells[match(colnames(counts), cells$barcode), ]

seurat_obj <- CreateSeuratObject(counts = counts,
                                 meta.data = cells,
                                 project = "Sample00632301")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")

spatial_mat <- as.matrix(cells[, c("x_centroid", "y_centroid")])
rownames(spatial_mat) <- cells$barcode

seurat_obj[["spatial"]] <- CreateDimReducObject(
  embeddings = spatial_mat,
  key = "spatial_",
  assay = DefaultAssay(seurat_obj)
)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "spatial",
  dims = 1:2,
  graph.name = "spatial_nn",
  k.param = 12
)
seurat_obj <- FindClusters(
  seurat_obj,
  graph.name = "spatial_nn",
  resolution = 0.4,
  algorithm = 1,
  name = "spatial_clusters"
)

SpatialDimPlot(seurat_obj, group.by = "seurat_clusters")
SpatialDimPlot(seurat_obj, group.by = "spatial_clusters")

saveRDS(seurat_obj, file = "Sample00632301_seurat.rds")
write.csv(seurat_obj@meta.data, file = "Sample00632301_metadata.csv")
