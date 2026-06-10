library(Seurat)

print("Loading Seurat object...")
seurat_obj <- readRDS("data/processed/spatial_with_timepoint.rds")
print("Seurat object loaded.")
print(seurat_obj)

print("Metadata columns:")
print(colnames(seurat_obj@meta.data))

print("Assays:")
print(names(seurat_obj@assays))

print("Reductions:")
print(names(seurat_obj@reductions))

print("First few rows of metadata:")
print(head(seurat_obj@meta.data))
