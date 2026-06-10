library(ggplot2)
library(Seurat)

binary_genes  <- c("Gja1", "Rhot1")
timepoints_to_plot <- c("Control", "3d", "14d", "28d")

target_celltypes <- c(
  "FAPs",
  "Lymphatic Endothelial Cells",
  "Endothelial Cells"
)

assay_use <- DefaultAssay(spatial)
feat_available <- rownames(GetAssayData(spatial, assay = assay_use, slot = "data"))
stopifnot(all(binary_genes %in% feat_available))

pdf_file <- "Spatial_AllCelltypes_Zooms_and_Binary.pdf"

build_ct_df <- function(spatial_obj, tp_label, ct_label) {
  cells_tp <- colnames(spatial_obj)[which(spatial_obj$timepoint == tp_label)]
  if (length(cells_tp) == 0) return(NULL)

  is_ct <- Idents(spatial_obj)[cells_tp] == ct_label
  cells_ct <- cells_tp[is_ct]
  if (length(cells_ct) == 0) return(NULL)

  df_ct <- data.frame(
    x_centroid = spatial_obj@meta.data[cells_ct, "x_centroid"],
    y_centroid = spatial_obj@meta.data[cells_ct, "y_centroid"],
    row.names  = cells_ct,
    check.names = FALSE
  )

  df_ct <- df_ct[
    is.finite(df_ct$x_centroid) & is.finite(df_ct$y_centroid),
    ,
    drop = FALSE
  ]

  df_ct$cell <- rownames(df_ct)
  return(df_ct)
}

binary_plot_zoom <- function(df_zoom, gene, x_range, y_range, tp, rep_name, ct_label){
  vals <- df_zoom[[gene]]
  df_zoom$binary <- ifelse(is.na(vals), FALSE, vals > 0)

  ggplot(df_zoom, aes(x = x_centroid, y = y_centroid, color = binary)) +
    geom_point(size = 0.32, alpha = 0.9) +
    coord_fixed(xlim = x_range, ylim = y_range) +
    scale_color_manual(
      values = c("TRUE" = "#60339b", "FALSE" = "#ffd014"),
      labels = c("FALSE" = paste0(gene,"-"), "TRUE" = paste0(gene,"+")),
      name = gene
    ) +
    labs(
      title = paste0(tp, " ", gene, "+/− ", ct_label, " ", rep_name),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}


pdf(pdf_file, width = 5, height = 5)

for (tp in timepoints_to_plot) {
  message("Timepoint: ", tp)

  for (ct_label in target_celltypes) {

    ct_df <- build_ct_df(spatial, tp, ct_label)
    if (is.null(ct_df) || nrow(ct_df) == 0) {
      warning("No cells of type ", ct_label, " for timepoint ", tp)
      next
    }

    expr_bin <- FetchData(spatial, vars = binary_genes, cells = ct_df$cell)
    expr_bin$cell <- rownames(expr_bin)

    df_ct_expr <- merge(ct_df, expr_bin, by = "cell", all.x = TRUE, sort = FALSE)

    for (rep_name in names(coords_list[[tp]])) {

      x_range <- coords_list[[tp]][[rep_name]]$x
      y_range <- coords_list[[tp]][[rep_name]]$y

      df_zoom <- subset(
        df_ct_expr,
        x_centroid >= x_range[1] & x_centroid <= x_range[2] &
          y_centroid >= y_range[1] & y_centroid <= y_range[2]
      )
      if (nrow(df_zoom) == 0) next

      p_base <- ggplot(df_zoom, aes(x = x_centroid, y = y_centroid)) +
        geom_point(size = 0.28, color = "#9179b8", alpha = 0.85) +
        coord_fixed(xlim = x_range, ylim = y_range) +
        labs(
          title = paste0(tp, " ", ct_label, " ", rep_name),
          x = NULL, y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)
        )

      print(p_base)

      for (gene in binary_genes) {
        p_bin <- binary_plot_zoom(df_zoom, gene, x_range, y_range, tp, rep_name, ct_label)
        print(p_bin)
      }
    }
  }
}

dev.off()

message("Finished PDF: ", pdf_file)
