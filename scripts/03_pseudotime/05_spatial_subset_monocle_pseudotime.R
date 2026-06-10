library(dplyr)
library(monocle3)
library(Seurat)
library(SeuratWrappers)

suppressMessages({
})


Idents(spatial) <- spatial@active.ident

target_clusters <- c(
  "FAPs",
  "Endothelial Cells",
  "Igf1+ Pro-Regenerative Fibroblasts",
  "Niche-Regulating Fibroblast",
  "Pericytes"
)

spatial_sub <- subset(spatial, idents = target_clusters)


cds <- as.cell_data_set(spatial_sub)

colData(cds)$seurat_cluster <- Idents(spatial_sub)

umap_mat <- Embeddings(spatial_sub, reduction = "umap")[, 1:2]

reducedDims(cds)$UMAP <- umap_mat

cds <- preprocess_cds(cds, num_dim = 50)

cds <- cluster_cells(cds, reduction_method = "UMAP")


cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)


fap_cells <- colnames(spatial_sub)[Idents(spatial_sub) == "FAPs"]

if (length(fap_cells) == 0) {
  stop("No cells with cluster == 'FAPs' found in spatial_sub.")
}

cds <- order_cells(cds, root_cells = fap_cells)


p_pseudo <- plot_cells(
  cds,
  reduction_method      = "UMAP",
  color_cells_by        = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups     = FALSE,
  label_leaves          = FALSE,
  label_branch_points   = FALSE,
  label_roots           = FALSE,
  graph_label_size      = 0   # no text bubbles / white circles
)

print(p_pseudo)

p_cluster <- plot_cells(
  cds,
  reduction_method      = "UMAP",
  color_cells_by        = "seurat_cluster",
  show_trajectory_graph = TRUE,
  label_cell_groups     = FALSE,
  label_leaves          = FALSE,
  label_branch_points   = FALSE,
  label_roots           = FALSE,
  graph_label_size      = 0
)

print(p_cluster)







de_spatial <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores = 4
)

sig_spatial <- de_spatial %>%
  dplyr::filter(q_value < 0.05) %>%
  dplyr::arrange(desc(morans_I))

sig_spatial$gene <- rownames(sig_spatial)
sig_spatial <- sig_spatial[, c("gene", setdiff(colnames(sig_spatial), "gene"))]


genes_to_plot <- head(sig_spatial$gene, 20)


get_expr_traj_spatial <- function(cds_branch, gene, seurat_obj, n_bins = 17, label = "Trajectory") {
  pseudo_vals <- monocle3::pseudotime(cds_branch)
  pseudo_vals <- pseudo_vals[!is.na(pseudo_vals)]

  if (length(pseudo_vals) == 0) {
    stop("No pseudotime values found in cds_branch.")
  }

  sub <- subset(seurat_obj, cells = names(pseudo_vals))

  sub$pseudo <- pseudo_vals[colnames(sub)]

  pt_min <- min(sub$pseudo, na.rm = TRUE)
  pt_max <- max(sub$pseudo, na.rm = TRUE)
  if (!is.finite(pt_min) || !is.finite(pt_max) || pt_min == pt_max) {
    stop("Pseudotime has zero or invalid range for this branch.")
  }

  breaks <- seq(pt_min, pt_max, length.out = n_bins + 1)

  sub$Pseudo_bin <- cut(
    sub$pseudo,
    breaks        = breaks,
    include.lowest = TRUE,
    labels        = FALSE
  )

  assay_use <- DefaultAssay(sub)  # "SCT" for you

  avg_exp <- Seurat::AverageExpression(
    sub,
    assay    = assay_use,
    features = gene,
    group.by = "Pseudo_bin"
  )

  mat <- avg_exp[[assay_use]]   # 1 x N matrix

  expr_vec <- as.numeric(mat[1, ])       # expression for bins in order
  bins_vec <- seq_along(expr_vec)        # 1..N

  df <- data.frame(
    Gene       = gene,
    Pseudo_bin = bins_vec,
    Expression = expr_vec,
    Trajectory = label,
    stringsAsFactors = FALSE
  )

  return(df)
}


plot_gene_traj_compare_onepdf <- function(
    cds_branch_endo,
    cds_branch_fibro,
    seurat_obj,
    genes,
    pdf_file = "~/Desktop/Spatial/pseudobin/endo_fibro/endo_vs_fibro_pseudotime.pdf",
    n_bins = 17
) {
  outdir <- dirname(pdf_file)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  pdf(pdf_file, width = 6, height = 4)
  on.exit(dev.off(), add = TRUE)

  for (gene in genes) {
    cat("Processing gene:", gene, "\n")

    if (!(gene %in% rownames(seurat_obj))) {
      cat("  -> Skipping:", gene, "(not found in Seurat)\n")
      next
    }
    if (!(gene %in% rownames(rowData(cds_branch_endo)))) {
      cat("  -> Skipping:", gene, "(not found in endothelial Monocle branch)\n")
      next
    }
    if (!(gene %in% rownames(rowData(cds_branch_fibro)))) {
      cat("  -> Skipping:", gene, "(not found in fibro Monocle branch)\n")
      next
    }

    df_endo <- get_expr_traj_spatial(
      cds_branch = cds_branch_endo,
      gene       = gene,
      seurat_obj = seurat_obj,
      n_bins     = n_bins,
      label      = "Endothelial"
    )

    if (!all(is.na(df_endo$Expression))) {
      p_endo <- ggplot(df_endo, aes(x = Pseudo_bin, y = Expression)) +
        geom_line() +
        geom_point() +
        theme_classic() +
        xlab("Pseudotime bin") +
        ylab(paste0(gene, " expression")) +
        ggtitle(paste0(gene, " – Endothelial lineage")) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text  = element_text(size = 10),
          axis.title = element_text(size = 12)
        )
      print(p_endo)   # page 1 for this gene
    } else {
      cat("  -> All NA expression for", gene, "in endothelial branch.\n")
    }

    df_fibro <- get_expr_traj_spatial(
      cds_branch = cds_branch_fibro,
      gene       = gene,
      seurat_obj = seurat_obj,
      n_bins     = n_bins,
      label      = "Fibrogenic"
    )

    if (!all(is.na(df_fibro$Expression))) {
      p_fibro <- ggplot(df_fibro, aes(x = Pseudo_bin, y = Expression)) +
        geom_line() +
        geom_point() +
        theme_classic() +
        xlab("Pseudotime bin") +
        ylab(paste0(gene, " expression")) +
        ggtitle(paste0(gene, " – Fibrogenic lineage")) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text  = element_text(size = 10),
          axis.title = element_text(size = 12)
        )
      print(p_fibro)  # page 2 for this gene
    } else {
      cat("  -> All NA expression for", gene, "in fibrogenic branch.\n")
    }
  }
}


genes_to_plot <- head(sig_spatial$gene, 20)

plot_gene_traj_spatial(
  cds_branch = cds,          # full trajectory on spatial_sub
  seurat_obj = spatial_sub,
  genes      = genes_to_plot,
  outdir     = "~/Desktop/Spatial/pseudobin/",
  n_bins     = 17
)

cds_endothelial <- choose_graph_segments(cds, clear_cds = FALSE)
cds_fibro <- choose_graph_segments(cds, clear_cds = FALSE)
endo_genes <- c("Fstl1", "Vegfa", "Angpt1", "Fgf2", "Gpihbp1", "Kdr")

fibro_genes <- c("Igf1", "Col1a1", "Col1a2", "Timp3", "Prg4", "Pi16")

plot_gene_traj_compare_onepdf(
  cds_branch_endo  = cds_endothelial,
  cds_branch_fibro = cds_fibro,
  seurat_obj       = spatial_sub,
  genes            = fibro_genes,
  pdf_file         = "~/Desktop/Spatial/pseudobin/endo_fibro/endo_vs_fibro_fibrogenes.pdf",
  n_bins           = 17
)
