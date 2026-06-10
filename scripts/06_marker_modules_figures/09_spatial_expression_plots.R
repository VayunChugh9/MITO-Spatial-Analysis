library(ggplot2)
library(Seurat)

stopifnot("timepoint" %in% colnames(spatial@meta.data))
stopifnot(all(c("x_centroid","y_centroid") %in% colnames(spatial@meta.data)))

target_level <- "28d"

cells_control <- colnames(spatial)[which(spatial$timepoint == target_level)]

idents_vec <- Idents(spatial)  # factor of active identities
idents_vec <- idents_vec[cells_control]

df <- data.frame(
  x_centroid  = spatial@meta.data[cells_control, "x_centroid"],
  y_centroid  = spatial@meta.data[cells_control, "y_centroid"],
  active_ident = idents_vec,
  sample_label = spatial@meta.data[cells_control, "sample_label"],
  row.names   = cells_control,
  check.names = FALSE
)

df <- df[is.finite(df$x_centroid) & is.finite(df$y_centroid), , drop = FALSE]

df$active_ident <- droplevels(df$active_ident)

p_control_idents <- ggplot(df, aes(x = x_centroid, y = y_centroid, color = active_ident)) +
  geom_point(size = 0.3) +
  coord_fixed() +
  labs(
    title = paste0(target_level, " Spatial Centroids"),
    x = "x_centroid", y = "y_centroid", color = "active.ident"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))  # << makes legend dots bigger


print(p_control_idents)




target_level <- "Control"
fap_label <- "FAPs"

cells_tp <- colnames(spatial)[which(spatial$timepoint == target_level)]

cells_faps <- cells_tp[Idents(spatial)[cells_tp] == fap_label]

df_faps <- data.frame(
  x_centroid  = spatial@meta.data[cells_faps, "x_centroid"],
  y_centroid  = spatial@meta.data[cells_faps, "y_centroid"],
  row.names   = cells_faps,
  check.names = FALSE
)
df_faps <- df_faps[is.finite(df_faps$x_centroid) & is.finite(df_faps$y_centroid), , drop = FALSE]

df_bg <- data.frame(
  x_centroid = spatial@meta.data[cells_tp, "x_centroid"],
  y_centroid = spatial@meta.data[cells_tp, "y_centroid"]
)
df_bg <- df_bg[is.finite(df_bg$x_centroid) & is.finite(df_bg$y_centroid), , drop = FALSE]


p_faps_only <- ggplot() +
  geom_point(data = df_bg, aes(x = x_centroid, y = y_centroid),
             size = 0.15, color = "grey85", alpha = 0.8) +
  geom_point(data = df_faps, aes(x = x_centroid, y = y_centroid),
             size = 0.35, color = "#9179b8") +
  coord_fixed() +
  labs(
    title = paste0(target_level, " ", fap_label),
    x = "x_centroid", y = "y_centroid"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
print(p_faps_only)








x_range <- c(5000, 8000)
y_range <- c(7000, 13000)

df_bg_zoom <- subset(df_bg,
                     x_centroid >= x_range[1] & x_centroid <= x_range[2] &
                       y_centroid >= y_range[1] & y_centroid <= y_range[2]
)
df_faps_zoom <- subset(df_faps,
                       x_centroid >= x_range[1] & x_centroid <= x_range[2] &
                         y_centroid >= y_range[1] & y_centroid <= y_range[2]
)

p_faps_zoom <- ggplot() +
  geom_point(
    data = df_bg_zoom,
    aes(x = x_centroid, y = y_centroid),
    size = 0.08, color = "grey85", alpha = 0.5
  ) +
  geom_jitter(
    data = df_faps_zoom,
    aes(x = x_centroid, y = y_centroid),
    width = 80, height = 80,  # adjust if you want less/more spread in the zoom
    size = 0.2,
    color = "#9179b8",
    alpha = 0.9
  ) +
  coord_fixed(xlim = x_range, ylim = y_range) +
  labs(
    title = paste0(target_level, " ", fap_label, " (x: ", x_range[1], "-", x_range[2],
                   ", y: ", y_range[1], "-", y_range[2], ")"),
    x = "x_centroid", y = "y_centroid"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )


print(p_faps_zoom)







gene <- "Kdr"
x_range <- c(4701, 8000)
y_range <- c(0, 7000)


target_level <- "28d"  # or whichever timepoint label you're working with

stopifnot(exists("spatial"))
stopifnot(all(c("x_centroid", "y_centroid", "timepoint") %in% colnames(spatial@meta.data)))

cells_tp <- colnames(spatial)[which(spatial$timepoint == target_level)]

df_tp <- data.frame(
  x_centroid = spatial@meta.data[cells_tp, "x_centroid"],
  y_centroid = spatial@meta.data[cells_tp, "y_centroid"],
  row.names  = cells_tp,
  check.names = FALSE
)
df_tp <- df_tp[is.finite(df_tp$x_centroid) & is.finite(df_tp$y_centroid), , drop = FALSE]
df_tp$cell <- rownames(df_tp)

assay_use <- DefaultAssay(spatial)
feat_available <- rownames(GetAssayData(spatial, assay = assay_use, slot = "data"))
if (!gene %in% feat_available) {
  stop(sprintf("Feature '%s' not found in assay '%s'.", gene, assay_use))
}

expr_df <- FetchData(spatial, vars = gene, cells = df_tp$cell)
expr_df$cell <- rownames(expr_df)

df_expr <- merge(df_tp, expr_df, by = "cell", all.x = TRUE, sort = FALSE)

df_expr_zoom <- subset(
  df_expr,
  x_centroid >= x_range[1] & x_centroid <= x_range[2] &
    y_centroid >= y_range[1] & y_centroid <= y_range[2]
)

p_myh2_zoom <- ggplot(df_expr_zoom, aes(x = x_centroid, y = y_centroid, color = !!sym(gene))) +
  geom_point(size = 0.3, alpha = 0.9, na.rm = TRUE) +
  coord_fixed(xlim = x_range, ylim = y_range) +
  scale_color_gradientn(
    colors = c("lightgrey", "red"),  # FeaturePlot-like palette
    name = gene
  ) +
  labs(
    title = paste0(target_level , " ", gene, " expression ",
                   "(x: ", x_range[1], "-", x_range[2], ", y: ", y_range[1], "-", y_range[2], ")"),
    x = "x_centroid", y = "y_centroid"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )


print(p_myh2_zoom)







gene <- "Rhot1"

stopifnot(exists("spatial"), exists("df_bg_zoom"), exists("df_faps_zoom"))
stopifnot(all(c("x_centroid","y_centroid") %in% colnames(df_faps_zoom)))

df_faps_zoom$cell <- rownames(df_faps_zoom)

assay_use <- DefaultAssay(spatial)
feat_available <- rownames(GetAssayData(spatial, assay = assay_use, slot = "data"))
if (!gene %in% feat_available) {
  stop(sprintf("Feature '%s' not found in assay '%s' (slot='data').", gene, assay_use))
}

cells_here <- df_faps_zoom$cell
expr_df <- FetchData(spatial, vars = gene, cells = cells_here)
expr_df$cell <- rownames(expr_df)

n_overlap <- sum(df_faps_zoom$cell %in% expr_df$cell)
message(sprintf("Matched %d of %d FAP zoom cells to Seurat cell names.", n_overlap, nrow(df_faps_zoom)))
if (n_overlap == 0) {
  warning("No IDs overlapped. Check slice/sample suffixes in metadata vs Seurat cell names.")
}

df_faps_zoom_expr <- merge(df_faps_zoom, expr_df, by = "cell", all.x = TRUE, sort = FALSE)

vals <- df_faps_zoom_expr[[gene]]
df_faps_zoom_expr$Rhot1_pos <- ifelse(is.na(vals), FALSE, vals > 0)

p_rhot1_binary_zoom <- ggplot() +

  geom_point(
    data = df_faps_zoom_expr,
    aes(x = x_centroid, y = y_centroid, color = Rhot1_pos),
    size = 0.25, alpha = 0.9, na.rm = TRUE
  ) +
  coord_fixed(xlim = x_range, ylim = y_range) +
  scale_color_manual(
    values = c("FALSE" = "#ffd014", "TRUE" = "#60339b"),
    labels  = c("Rhot1-", "Rhot1+"),
    name    = gene
  ) +
  labs(
    title = paste0(target_level, " FAPs ", gene, "+ vs ", gene, "- (x: ",
                   x_range[1], "-", x_range[2], ", y: ", y_range[1], "-", y_range[2], ")"),
    x = "x_centroid",
    y = "y_centroid"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    plot.title  = element_text(hjust = 0.5)
  )


print(p_rhot1_binary_zoom)
