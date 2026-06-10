library(ggplot2)
library(Seurat)

stopifnot(exists("spatial"))
stopifnot("timepoint" %in% colnames(spatial@meta.data))
stopifnot(all(c("x_centroid","y_centroid") %in% colnames(spatial@meta.data)))

celltype_col <- NULL   # e.g., "celltype"  (leave NULL to use active.ident)

build_df_timepoint <- function(obj, tp, celltype_col = NULL) {
  cells_tp <- colnames(obj)[which(obj$timepoint == tp)]

  if (length(cells_tp) == 0) stop(paste0("No cells found for timepoint: ", tp))

  if (!is.null(celltype_col)) {
    stopifnot(celltype_col %in% colnames(obj@meta.data))
    ct <- obj@meta.data[cells_tp, celltype_col]
  } else {
    ct <- Idents(obj)[cells_tp]   # active.ident
  }

  df <- data.frame(
    x_centroid = obj@meta.data[cells_tp, "x_centroid"],
    y_centroid = obj@meta.data[cells_tp, "y_centroid"],
    cell_type  = ct,
    row.names  = cells_tp,
    check.names = FALSE
  )

  df <- df[is.finite(df$x_centroid) & is.finite(df$y_centroid), , drop = FALSE]
  df$cell_type <- droplevels(factor(df$cell_type))
  df
}

plot_timepoint_xy <- function(df, tp, flip_y = FALSE, pt_size = 0.25) {
  p <- ggplot(df, aes(x = x_centroid, y = y_centroid, color = cell_type)) +
    geom_point(size = pt_size) +
    coord_fixed() +
    labs(
      title = paste0(tp, " Spatial Celltypes"),
      x = "x_centroid", y = "y_centroid", color = "cell type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  if (flip_y) p <- p + scale_y_reverse()
  p
}

tps <- levels(factor(spatial$timepoint))  # should give Control, 3d, 14d, 28d (order depends on factor levels)
plots <- setNames(vector("list", length(tps)), tps)

for (tp in tps) {
  df_tp <- build_df_timepoint(spatial, tp, celltype_col = celltype_col)
  plots[[tp]] <- plot_timepoint_xy(df_tp, tp, flip_y = FALSE, pt_size = 0.25)
  print(plots[[tp]])
}
