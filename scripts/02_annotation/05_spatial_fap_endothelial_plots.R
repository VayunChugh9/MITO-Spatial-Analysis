library(ggplot2)
library(Seurat)

pathway_genes <- c(
  "Pdgfra",                        # FAP marker
  "Bnip3", "Bnip3l", "Higd1b",     # Hypoxia/Metabolic stress proxies
  "Vegfa", "Igf1", "Fgf2",         # PI3K/mTOR & Angiogenesis drivers
  "Nos2",                          # NO synthesis proxy
  "Kdr", "Eng", "Tie1", "Sox17",   # Endothelial differentiation
  "Vwf", "Gja1"                    # Mature Endothelial markers
)


spatial_pathway_genes <- c(
  "Sox17", "Tie1", "Kdr", "Eng",      # early / angiogenic / EC commitment
  "Vegfa", "Nos2", "Gja1",            # functional/activation & coupling
  "Vwf",                               # mature EC marker
  "Pdgfra", "Igf1", "Fgf2",            # non-EC / paracrine (kept later)
  "Bnip3", "Bnip3l", "Higd1b"          # hypoxia/stress (kept last)
)

spatial$timepoint <- factor(spatial$timepoint,
                            levels = c("Control", "3d", "14d", "28d"))

faps_spatial <- subset(spatial, idents = "FAPs")
endo_spatial <- subset(spatial, idents = "Endothelial Cells")

DotPlot(faps_spatial,
        features = spatial_pathway_genes,
        group.by = "timepoint") +
  ggtitle("FAP to EC Pathway Expression in Spatial FAPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

DotPlot(endo_spatial,
        features = spatial_pathway_genes,
        group.by = "timepoint") +
  ggtitle("FAP to EC Pathway Expression in Spatial Endothelial Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())



spatial_pathway_genes <- c("Pdgfra", "Bnip3", "Bnip3l", "Higd1b", "Vegfa",
                           "Igf1", "Fgf2", "Nos2", "Kdr", "Eng", "Tie1",
                           "Sox17", "Vwf", "Gja1")

spatial$timepoint <- factor(spatial$timepoint,
                            levels = c("Control", "3d", "14d", "28d"))

timepoints <- levels(spatial$timepoint)

pdf_path <- "data/processed/Spatial_Fap_Endo_UMAPs.pdf"
pdf(pdf_path, width = 8, height = 6)

cat("Generating PDF... this will take a moment as it plots all cells.\n")

for (gene in spatial_pathway_genes) {

  if (gene %in% rownames(spatial)) {

    for (tp in timepoints) {
      tp_subset <- subset(spatial, timepoint == tp)

      p <- FeaturePlot(tp_subset, features = gene, order = TRUE, pt.size = 0.5) +
        ggtitle(paste(gene, "Expression at", tp)) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

      print(p)
    }
  } else {
    cat(paste("Warning: Gene", gene, "not found in dataset. Skipping...\n"))
  }
}

dev.off()
cat(paste("Done! PDF saved to:", pdf_path, "\n"))



spatial$timepoint <- factor(spatial$timepoint,
                            levels = c("Control", "3d", "14d", "28d"))

pdf_path <- "data/processed/Spatial_DimPlots_by_Timepoint.pdf"
pdf(pdf_path, width = 8, height = 6)

timepoints <- levels(spatial$timepoint)

cat("Generating Cell Type UMAP PDF...\n")

for (tp in timepoints) {
  tp_subset <- subset(spatial, timepoint == tp)

  p <- DimPlot(tp_subset,

               repel = TRUE,
               label.size = 3,
               pt.size = 0.5) +
    ggtitle(paste(tp)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          legend.position = "right")

  print(p)
}

dev.off()
cat(paste("Done! DimPlot PDF saved to:", pdf_path, "\n"))



faps_subset$timepoint <- factor(faps_subset$timepoint, levels = c("Control", "3d", "14d", "28d"))

bnip3l_counts <- GetAssayData(faps_subset, assay = DefaultAssay(faps_subset), layer = "data")["Bnip3l", ]

faps_subset$Bnip3l_Status <- ifelse(bnip3l_counts > 0, "Bnip3l_Pos", "Bnip3l_Neg")

faps_subset$Time_and_Bnip3l <- paste(faps_subset$timepoint, faps_subset$Bnip3l_Status, sep = " - ")

faps_subset$Time_and_Bnip3l <- factor(faps_subset$Time_and_Bnip3l, levels = c(
  "Control - Bnip3l_Neg", "Control - Bnip3l_Pos",
  "3d - Bnip3l_Neg",      "3d - Bnip3l_Pos",
  "14d - Bnip3l_Neg",     "14d - Bnip3l_Pos",
  "28d - Bnip3l_Neg",     "28d - Bnip3l_Pos"
))

DotPlot(faps_subset, features = c("Gja1", "Vwf", "Kdr"), group.by = "Time_and_Bnip3l", dot.scale = 8) +
  ggtitle("EC Fate Markers in Bnip3l+ vs Bnip3l- FAPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylab("Timepoint & Hypoxia (Bnip3l) Status") +
  xlab("Gene")
