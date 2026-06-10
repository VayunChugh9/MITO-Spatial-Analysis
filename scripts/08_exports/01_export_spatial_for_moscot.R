library(data.table)
library(Matrix)
library(Seurat)

print("Loading Seurat object...")
seurat_obj <- readRDS("data/processed/spatial_with_timepoint.rds")

meta <- seurat_obj@meta.data
meta$barcode <- rownames(meta)

if (!all(c("x_centroid", "y_centroid") %in% colnames(meta))) {
    stop("x_centroid or y_centroid missing in metadata")
}

print("Exporting metadata...")
write.csv(meta, "spatial_metadata.csv", row.names = FALSE)

print("Exporting counts...")
counts <- GetAssayData(seurat_obj, slot = "counts")
writeMM(counts, "spatial_counts.mtx")
write.table(rownames(counts), "spatial_features.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(counts), "spatial_barcodes.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)

print("Export complete.")
