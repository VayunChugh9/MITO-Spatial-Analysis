library(Seurat)
library(SeuratDisk)

base_dir <- "data/processed/Spatiotemporal_Analysis/CyteType Temporal Problem"
results_dir <- file.path(base_dir, "results")

rds_path <- "data/processed/spatial_cyte_type.rds"

h5ad_path <- file.path(results_dir, "spatial_cyte_type.h5ad")

if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

print(paste(rep("=", 60), collapse = ""))
print("Step 0: Converting RDS to h5ad")
print(paste("Input:", rds_path))
print(paste("Output:", h5ad_path))
print(paste(rep("=", 60), collapse = ""))

if (!file.exists(rds_path)) stop(paste("RDS not found at:", rds_path))

cat("\nLoading RDS...\n")
spatial <- readRDS(rds_path)

cat("Processing Metadata (CyteType_Clusters -> celltype)...\n")
if ("CyteType_Clusters" %in% colnames(spatial@meta.data)) {
    spatial$celltype <- as.character(spatial@meta.data$CyteType_Clusters)
    Idents(spatial) <- spatial$celltype
    cat(paste("  - Mapped", length(unique(spatial$celltype)), "clusters to 'celltype'\n"))
} else {
    stop("ERROR: Metadata column 'CyteType_Clusters' not found!")
}

spatial@images <- list()

cat("Converting to h5ad...\n")
temp_h5seurat <- file.path(results_dir, "temp_conversion.h5seurat")

tryCatch({
    SaveH5Seurat(spatial, filename = temp_h5seurat, overwrite = TRUE, verbose = FALSE)
    Convert(temp_h5seurat, dest = "h5ad", overwrite = TRUE)

    generated_file <- sub(".h5seurat", ".h5ad", temp_h5seurat)

    if (file.exists(generated_file)) {
        file.rename(generated_file, h5ad_path)
        cat(paste("  - Success: Saved to", h5ad_path, "\n"))
    } else {
        stop("Conversion generated no output file.")
    }
}, finally = {
    if (file.exists(temp_h5seurat)) file.remove(temp_h5seurat)
})
