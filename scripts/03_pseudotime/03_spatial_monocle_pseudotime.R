library(dplyr)
library(FNN)
library(ggplot2)
library(monocle3)
library(paletteer)
library(Seurat)
library(SeuratWrappers)

Idents(spatial) <- spatial@active.ident
clusters <- Idents(spatial)

set.seed(123)  # for reproducibility

target_total <- 40000   # adjust up/down based on your RAM

tab <- table(clusters)
cluster_prop <- tab / sum(tab)

target_per_cluster <- ceiling(cluster_prop * target_total)
target_per_cluster <- pmin(target_per_cluster, tab)

cells_by_cluster <- split(names(clusters), clusters)

cells_sub <- unlist(
  lapply(names(cells_by_cluster), function(cl){
    pool <- cells_by_cluster[[cl]]
    n_target <- target_per_cluster[[cl]]
    if (length(pool) <= n_target) {
      pool
    } else {
      sample(pool, n_target)
    }
  }),
  use.names = FALSE
)

length(cells_sub)  # how many cells you actually kept

spatial_sub <- subset(spatial, cells = cells_sub)

if (!"FAPs" %in% levels(Idents(spatial_sub))) {
  stop("FAPs cluster not present after downsampling – increase target_total or sample FAPs separately.")
}

suppressPackageStartupMessages({
})

cds_sub <- as.cell_data_set(spatial_sub)

colData(cds_sub)$seurat_cluster <- Idents(spatial_sub)

umap_sub <- Embeddings(spatial_sub, reduction = "umap")
umap_sub <- umap_sub[colnames(cds_sub), ]   # align

reducedDims(cds_sub)$UMAP <- umap_sub

cds_sub <- preprocess_cds(cds_sub, num_dim = 30)
cds_sub <- cluster_cells(cds_sub, reduction_method = "UMAP")

cds_sub <- learn_graph(
  cds_sub,
  use_partition = TRUE,
  close_loop = FALSE
)

fap_cells_sub <- colnames(spatial_sub)[Idents(spatial_sub) == "FAPs"]
if (length(fap_cells_sub) == 0) {
  stop("No FAPs in spatial_sub – check sampling.")
}

cds_sub <- order_cells(
  cds_sub,
  reduction_method = "UMAP",
  root_cells = fap_cells_sub
)

pt_sub <- pseudotime(cds_sub)
head(pt_sub)


pca_full <- Embeddings(spatial, reduction = "pca")
pca_sub  <- pca_full[colnames(spatial_sub), ]

pt_sub <- pt_sub[colnames(spatial_sub)]

k <- 10  # neighbors

knn_res <- get.knnx(
  data  = pca_sub,
  query = pca_full,
  k     = k
)

pt_sub_vec <- as.numeric(pt_sub)
names(pt_sub_vec) <- colnames(spatial_sub)

pt_full <- sapply(1:nrow(knn_res$nn.index), function(i){
  idx <- knn_res$nn.index[i, ]
  neighbor_cells <- rownames(pca_sub)[idx]
  mean(pt_sub_vec[neighbor_cells], na.rm = TRUE)
})

names(pt_full) <- rownames(pca_full)

spatial$pseudotime <- pt_full[colnames(spatial)]

summary(spatial$pseudotime)

FeaturePlot(spatial, features = "pseudotime", reduction = "umap")


p <- plot_cells(
  cds_sub,
  color_cells_by        = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups     = FALSE,
  label_leaves          = FALSE,
  label_branch_points   = FALSE,
  graph_label_size      = 0
)

p + scale_color_gradientn(
  name    = "pseudotime",   # legend title
  colours = as.character(
    paletteer_c("ggthemes::Orange-Blue Diverging", 100)
  )
)






plot_cells(
  cds_sub,
  reduction_method = "UMAP",
  color_cells_by = "seurat_cluster",   # we stored this earlier
  show_trajectory_graph = TRUE,
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  label_roots = TRUE
)
