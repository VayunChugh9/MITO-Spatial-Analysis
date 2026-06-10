library(dplyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)

Idents(spatial) <- spatial@active.ident
clusters <- Idents(spatial)                   # factor of cluster labels
umap_mat <- Embeddings(spatial, reduction = "umap")

DefaultAssay(spatial) <- "SCT"                # make sure SCT exists

sce <- SingleCellExperiment(
  assays = list(
    logcounts = spatial[["SCT"]]@data        # log-normalized SCT values
  ),
  reducedDims = list(
    UMAP = umap_mat
  )
)

root_cluster <- "FAPs"

if (!root_cluster %in% levels(clusters)) {
  stop(paste("Root cluster", root_cluster, "not found in Idents(spatial).",
             "Available:", paste(levels(clusters), collapse = ", ")))
}

colLabels(sce) <- clusters

sce <- slingshot(
  sce,
  clusterLabels = colLabels(sce),
  reducedDim    = "pca",
  start.clus    = root_cluster   # <---- key line
)

lin_ind <- 1
pt_mat  <- slingPseudotime(sce, na = FALSE)
pt      <- pt_mat[, lin_ind]

spatial[[paste0("pseudotime_lineage", lin_ind)]] <- pt

df_pts <- data.frame(
  UMAP1      = umap_mat[, 1],
  UMAP2      = umap_mat[, 2],
  pseudotime = pt
)

crv <- slingCurves(sce)[[lin_ind]]

df_curve <- data.frame(
  UMAP1 = crv$s[, 1],
  UMAP2 = crv$s[, 2],
  ord   = seq_len(nrow(crv$s))
)

ggplot(df_pts, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = pseudotime), size = 0.5, na.rm = TRUE) +
  geom_path(
    data  = df_curve[order(df_curve$ord), ],
    aes(UMAP1, UMAP2),
    inherit.aes = FALSE,
    linewidth = 1,
    color = "black"
  ) +
  scale_color_viridis_c(name = "Pseudotime") +
  coord_equal() +
  theme_classic() +
  labs(
    title = paste0("Slingshot pseudotime (Lineage ", lin_ind,
                   "), root = ", root_cluster),
    x = "UMAP 1",
    y = "UMAP 2"
  )



for (lin_ind in seq_len(n_lin)) {
  spatial[[paste0("pseudotime_lineage", lin_ind)]] <- pt_mat[, lin_ind]
}

head(spatial@meta.data[, grepl("^pseudotime_lineage", colnames(spatial@meta.data))])






library(grid)        # for arrow()

umap_mat <- Embeddings(spatial, "umap")

df_all <- data.frame(
  UMAP1   = umap_mat[, 1],
  UMAP2   = umap_mat[, 2],
  cluster = Idents(spatial)
)

df_cent <- df_all %>%
  group_by(cluster) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  )

lin_list <- slingLineages(sce)
n_lin    <- length(lin_list)

pdf("Lineage_MST_UMAPs.pdf", width = 8, height = 6)

for (i in seq_len(n_lin)) {

  lin_clusters <- lin_list[[i]]           # vector of cluster names
  lin_label    <- if (!is.null(names(lin_list))) {
    names(lin_list)[i]
  } else {
    paste0("Lineage ", i)
  }

  if (length(lin_clusters) > 1) {
    df_edges <- do.call(rbind, lapply(seq_len(length(lin_clusters) - 1), function(k) {
      data.frame(from = lin_clusters[k], to = lin_clusters[k + 1])
    }))

    df_edges <- df_edges %>%
      left_join(df_cent, by = c("from" = "cluster")) %>%
      rename(x1 = UMAP1, y1 = UMAP2) %>%
      left_join(df_cent, by = c("to" = "cluster")) %>%
      rename(x2 = UMAP1, y2 = UMAP2)
  } else {
    df_edges <- data.frame(
      from = character(0), to = character(0),
      x1 = numeric(0), y1 = numeric(0),
      x2 = numeric(0), y2 = numeric(0)
    )
  }

  p <- ggplot() +
    geom_point(
      data  = df_all,
      aes(UMAP1, UMAP2, color = cluster),
      size  = 0.4,
      alpha = 0.6
    ) +

    geom_point(
      data  = df_cent,
      aes(UMAP1, UMAP2),
      size  = 2.5,
      color = "black"
    ) +

    geom_segment(
      data  = df_edges,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color    = "black",
      linewidth = 1.2,
      arrow    = arrow(length = unit(0.24, "cm"), type = "closed")
    ) +

    geom_text_repel(
      data = df_cent,
      aes(UMAP1, UMAP2, label = cluster),
      size         = 3,
      max.overlaps = 100
    ) +

    coord_equal() +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title      = element_text(hjust = 0.5)  # center the title
    ) +
    labs(
      title = paste0("MST on UMAP - ", lin_label),
      x = "UMAP1",
      y = "UMAP2",
      color = "Cluster"
    )

  print(p)  # one page in the PDF
}

dev.off()
