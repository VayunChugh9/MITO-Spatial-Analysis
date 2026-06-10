library(data.table)
library(ggplot2)
library(Seurat)

spatial <- readRDS("data/processed/spatial_with_timepoint.rds")

suppressPackageStartupMessages({
})

files <- data.frame(
  suffix = c(1, 2, 3, 4),
  path   = c(
    "data/processed/Sample00632261/cells61.csv",
    "data/processed/Sample00632262/cells62.csv",
    "data/processed/Sample00632301/cells01.csv",
    "data/processed/Sample00632302/cells02.csv"
  ),
  stringsAsFactors = FALSE
)

files$sample_label <- c("00632261 (Control)", "00632262", "00632301", "00632302")

read_one <- function(path, suffix, sample_label) {
  dt <- fread(path)
  stopifnot(all(c("cell_id","x_centroid","y_centroid") %in% names(dt)))
  dt[, cell_id := as.character(cell_id)]
  dt[, x_centroid := as.numeric(x_centroid)]
  dt[, y_centroid := as.numeric(y_centroid)]
  dt[, cell_id_with_suffix := paste0(cell_id, "_", suffix)]
  dt[, sample_suffix := suffix]
  dt[, sample_file   := basename(path)]
  dt[, sample_label  := sample_label]
  dt[]
}

meta_all <- rbindlist(Map(read_one, files$path, files$suffix, files$sample_label), fill = TRUE)

meta_all <- meta_all[!duplicated(cell_id_with_suffix)]
rn <- meta_all$cell_id_with_suffix
meta_all_df <- as.data.frame(meta_all)
rownames(meta_all_df) <- rn

stopifnot(exists("spatial"), inherits(spatial, "Seurat"))
cells_in_obj <- colnames(spatial)

drop_cols <- c("cell_id", "cell_id_with_suffix")  # we already use these as rownames/keys
keep_cols <- setdiff(colnames(meta_all_df), drop_cols)
meta_to_add <- meta_all_df[cells_in_obj, keep_cols, drop = FALSE]

spatial <- AddMetaData(spatial, metadata = meta_to_add)

n_added <- sum(c("x_centroid","y_centroid") %in% colnames(spatial@meta.data))
message("Added centroid columns? ", n_added == 2)

color_by <- if ("seurat_clusters" %in% colnames(spatial@meta.data)) "seurat_clusters" else "sample_label"

df <- FetchData(spatial, vars = c("x_centroid","y_centroid", color_by))
df <- df[is.finite(df$x_centroid) & is.finite(df$y_centroid), , drop = FALSE]

p <- ggplot(df, aes(x = x_centroid, y = y_centroid, color = .data[[color_by]])) +
  geom_point(size = 0.25) +
  coord_fixed() +
  labs(x = "x_centroid", y = "y_centroid", color = color_by, title = "Spatial centroids") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))


print(p)
