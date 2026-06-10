library(Seurat)
library(SeuratDisk)

rds_path <- "data/processed/spatial_with_timepoint.rds"
h5ad_path <- "data/processed/spatial_with_timepoint.h5ad"

if (!file.exists(rds_path)) {
    stop(paste("RDS file not found at:", rds_path))
}

cat(paste(rep("=", 60), collapse=""), "\n")
cat("Step 0: Converting RDS to h5ad (Fixed & Optimized)\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("\nLoading Seurat object from RDS...\n")
spatial <- readRDS(rds_path)

cat("  - Seurat object loaded successfully\n")
cat("  - Number of cells:", ncol(spatial), "\n")
cat("  - Number of features:", nrow(spatial), "\n")

cat("\nChecking Cell Type Annotations...\n")

active_levels <- levels(spatial)
n_active <- length(active_levels)
cat(sprintf("  - Active Ident levels found: %d\n", n_active))

cat("  - \033[33mACTION:\033[0m Overwriting 'celltype' with Active Ident labels (as text).\n")
spatial$celltype <- as.character(Idents(spatial))

n_meta <- length(unique(spatial@meta.data$celltype))
if (n_meta == n_active) {
     cat("  - \033[32mSUCCESS:\033[0m Metadata 'celltype' now contains", n_meta, "text labels.\n")
} else {
     warning("Mismatch remaining between active ident and metadata.")
}

cat("\nCleaning Spatial Object...\n")
if (length(spatial@images) > 0) {
    cat("  - Found spatial images:", paste(names(spatial@images), collapse=", "), "\n")
    cat("  - \033[33mACTION:\033[0m Removing all spatial FOVs to prevent export issues.\n")

    spatial@images <- list()

    cat("  - Spatial images removed. Object is now lighter for conversion.\n")
} else {
    cat("  - No spatial images found to remove.\n")
}

if ("timepoint" %in% colnames(spatial@meta.data)) {
    timepoints <- unique(spatial@meta.data$timepoint)
    cat("  - Timepoints found:", paste(timepoints, collapse=", "), "\n")
} else {
    warning("'timepoint' column not found in metadata")
}

cat("\nConverting to h5ad format using SeuratDisk...\n")

output_dir <- dirname(h5ad_path)

clean_name <- sub("\\.[^.]*$", "", basename(h5ad_path))
temp_h5seurat <- file.path(output_dir, paste0("temp_", clean_name, ".h5seurat"))

if (file.exists(temp_h5seurat)) {
    file.remove(temp_h5seurat)
}

tryCatch({
    cat("  - Step 1: Saving as h5Seurat format...\n")
    SaveH5Seurat(spatial, filename = temp_h5seurat, overwrite = TRUE, verbose = FALSE)

    if (!file.exists(temp_h5seurat)) {
        stop("Failed to create h5Seurat file")
    }
    cat("  - ✓ h5Seurat file created successfully\n")

    cat("  - Step 2: Converting to h5ad format...\n")

    conversion_success <- FALSE
    tryCatch({
        suppressWarnings({
            Convert(temp_h5seurat, dest = "h5ad", overwrite = TRUE)
        })
        conversion_success <- TRUE
        cat("  - Convert() completed without errors\n")
    }, error = function(e) {
        cat("  - ❌ Convert() threw an error:\n")
        cat("     ", conditionMessage(e), "\n")
        conversion_success <- FALSE
    })

    if (!conversion_success) {
        stop("Convert() function failed")
    }

    expected_file_1 <- sub("\\.h5seurat$", ".h5ad", temp_h5seurat)
    expected_file_2 <- paste0(temp_h5seurat, ".h5ad")

    converted_file <- NULL
    if (file.exists(expected_file_1)) {
        converted_file <- expected_file_1
    } else if (file.exists(expected_file_2)) {
        converted_file <- expected_file_2
    } else {
        cat("  - Standard names not found, scanning directory...\n")
        pattern <- paste0("^temp_", clean_name, ".*\\.h5ad$")
        found <- list.files(output_dir, pattern = pattern, full.names = TRUE)
        if(length(found) > 0) converted_file <- found[1]
    }

    if (is.null(converted_file)) {
        cat("  - Files in output directory:\n")
        print(list.files(output_dir, pattern = "temp_", full.names = TRUE))
        stop("Conversion failed - h5ad file not found.")
    }

    cat("  - ✓ h5ad file created:", converted_file, "\n")

    if (file.exists(h5ad_path)) {
        file.remove(h5ad_path)
    }
    file.rename(converted_file, h5ad_path)

    if (!file.exists(h5ad_path)) {
        stop("Failed to move converted file to final location")
    }

    cat("  - ✓ File saved to final location:", h5ad_path, "\n")

}, error = function(e) {
    cat("\n❌ Error during conversion:\n")
    cat("  ", conditionMessage(e), "\n")
    stop(conditionMessage(e))
}, finally = {
    if (file.exists(temp_h5seurat)) {
        file.remove(temp_h5seurat)
    }
})

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("Conversion complete!\n")
cat("Output saved to:", h5ad_path, "\n")
cat(paste(rep("=", 60), collapse=""), "\n")
