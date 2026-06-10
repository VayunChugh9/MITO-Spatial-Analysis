library(ggplot2)
library(Seurat)

gene <- "Gja1"
timepoints_to_plot <- c("Control", "3d", "14d", "28d")

target_celltypes <- c(
  "Activated Satellite Cells",
  "Quiescent Satellite Cells",
  "Myod1+/Myog+ Differentiating Myoblasts"
)

sub <- subset(spatial, celltype %in% target_celltypes)
sub$timepoint <- factor(sub$timepoint, levels = timepoints_to_plot)

pdf("Gja1_SatCells_Myoblasts_FeaturePlot.pdf", width = 6, height = 5)

for (ct_label in target_celltypes) {
  ct_sub <- subset(sub, celltype == ct_label)

  for (tp in timepoints_to_plot) {
    tp_sub <- subset(ct_sub, timepoint == tp)
    if (ncol(tp_sub) == 0) next

    p <- FeaturePlot(tp_sub, features = gene, order = TRUE) +
      ggtitle(paste0(tp, " - ", ct_label)) +
      theme(plot.title = element_text(hjust = 0.5))

    print(p)
  }
}

dev.off()
message("Done: Gja1_SatCells_Myoblasts_FeaturePlot.pdf")
