library(ggplot2)
library(Seurat)

library(rlang)   # for sym()


genes_to_plot <- c("Gja1", "Fstl1", "Vegfa", "Angpt1", "Fgf2", "Gpihbp1", "Kdr", 'Myh2', 'Myh3', 'Myh7')

timepoints_to_plot <- c("Control", "3d", "14d", "28d")

pdf_file <- "Spatial_Feature_Zooms_AllGenes_AllTimepoints.pdf"


coords_list <- list(
  "Control" = list(
    "Replicate 1" = list(
      x = c(0, 4999),
      y = c(1000, 6000)
    ),
    "Replicate 2" = list(
      x = c(5000, 8000),
      y = c(7000, 13000)
    ),
    "Replicate 3" = list(
      x = c(5500, 10000),
      y = c(0, 5000)
    )
  ),
  "3d" = list(
    "Replicate 1" = list(
      x = c(2500, 6000),
      y = c(5000, 8000)
    ),
    "Replicate 2" = list(
      x = c(0, 2500),
      y = c(0, 5000)
    ),
    "Replicate 3" = list(
      x = c(3000, 8000),
      y = c(1000, 5000)
    )
  ),
  "14d" = list(
    "Replicate 1" = list(
      x = c(0, 4000),
      y = c(2000, 6000)
    ),
    "Replicate 2" = list(
      x = c(4001, 6500),
      y = c(4100, 8000)
    ),
    "Replicate 3" = list(
      x = c(4300, 8000),
      y = c(0, 3000)
    )
  ),
  "28d" = list(
    "Replicate 1" = list(
      x = c(500, 4700),
      y = c(0, 6000)
    ),
    "Replicate 2" = list(
      x = c(0, 6000),
      y = c(6000, 12000)
    ),
    "Replicate 3" = list(
      x = c(4701, 8000),
      y = c(0, 7000)
    )
  )
)


stopifnot(exists("spatial"))
stopifnot(all(c("x_centroid", "y_centroid", "timepoint") %in% colnames(spatial@meta.data)))

assay_use <- DefaultAssay(spatial)
feat_available <- rownames(GetAssayData(spatial, assay = assay_use, slot = "data"))


build_tp_df <- function(spatial_obj, tp_label) {
  cells_tp <- colnames(spatial_obj)[which(spatial_obj$timepoint == tp_label)]
  df_tp <- data.frame(
    x_centroid = spatial_obj@meta.data[cells_tp, "x_centroid"],
    y_centroid = spatial_obj@meta.data[cells_tp, "y_centroid"],
    row.names  = cells_tp,
    check.names = FALSE
  )
  df_tp <- df_tp[
    is.finite(df_tp$x_centroid) & is.finite(df_tp$y_centroid),
    ,
    drop = FALSE
  ]
  df_tp$cell <- rownames(df_tp)
  df_tp
}


pdf(pdf_file, width = 5, height = 5)

for (tp in timepoints_to_plot) {
  if (!tp %in% names(coords_list)) {
    warning(sprintf("No coordinates defined for timepoint '%s'; skipping.", tp))
    next
  }

  message("Processing timepoint: ", tp)

  df_tp <- build_tp_df(spatial, tp)

  if (nrow(df_tp) == 0) {
    warning(sprintf("No cells found for timepoint '%s'; skipping.", tp))
    next
  }

  for (gene in genes_to_plot) {
    if (!gene %in% feat_available) {
      warning(sprintf("Feature '%s' not found in assay '%s'; skipping.", gene, assay_use))
      next
    }

    message("  Gene: ", gene)

    expr_df <- FetchData(spatial, vars = gene, cells = df_tp$cell)
    expr_df$cell <- rownames(expr_df)

    df_expr <- merge(df_tp, expr_df, by = "cell", all.x = TRUE, sort = FALSE)

    gene_min <- min(df_expr[[gene]], na.rm = TRUE)
    gene_max <- max(df_expr[[gene]], na.rm = TRUE)

    if (!is.finite(gene_min) || !is.finite(gene_max) || gene_min == gene_max) {
      gene_min <- 0
      gene_max <- max(df_expr[[gene]], na.rm = TRUE)
      if (!is.finite(gene_max) || gene_max == 0) {
        gene_max <- 1
      }
    }

    for (rep_name in c("Replicate 1", "Replicate 2", "Replicate 3")) {
      if (!rep_name %in% names(coords_list[[tp]])) {
        warning(sprintf("No coordinates for %s at timepoint '%s'; skipping.", rep_name, tp))
        next
      }

      ranges <- coords_list[[tp]][[rep_name]]
      x_range <- ranges$x
      y_range <- ranges$y

      df_expr_zoom <- subset(
        df_expr,
        x_centroid >= x_range[1] & x_centroid <= x_range[2] &
          y_centroid >= y_range[1] & y_centroid <= y_range[2]
      )

      if (nrow(df_expr_zoom) == 0) {
        warning(sprintf(
          "No cells in zoom window for %s, %s, %s; skipping plot.",
          tp, gene, rep_name
        ))
        next
      }

      p_zoom <- ggplot(df_expr_zoom, aes(x = x_centroid, y = y_centroid, color = !!sym(gene))) +
        geom_point(size = 0.3, alpha = 0.9, na.rm = TRUE) +
        coord_fixed(xlim = x_range, ylim = y_range) +
        scale_color_gradientn(
          colors = c("lightgrey", "red"),
          limits = c(gene_min, gene_max),
          name = gene
        ) +
        labs(
          title = paste0(
            tp, " ", gene, " ", rep_name, " (x: ",
            x_range[1], "-", x_range[2], ", y: ",
            y_range[1], "-", y_range[2], ")"
          ),
          x = "x_centroid",
          y = "y_centroid"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )

      print(p_zoom)
    }
  }
}

dev.off()

message("Finished writing PDF: ", pdf_file)
