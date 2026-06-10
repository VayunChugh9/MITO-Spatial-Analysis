library(dplyr)
library(ggplot2)
library(tidyr)

spatial <- readRDS('data/processed/spatial_with_timepoint.rds')
spatial_cyte_type <- readRDS("data/processed/spatial_cyte_type.rds")


FeaturePlot(spatial, features = "Gja1", reduction = "umap", order = TRUE)


VlnPlot(
  spatial,
  features = "Gja1",
  group.by = "timepoint",
  pt.size = 0,  # removes points
  cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")  # pick colors you like
) +

  ggtitle("Gja1 Expression by Timepoint") +
  ylab("Expression Level") +
  xlab("Timepoint") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # centered title
    axis.text.x = element_text(angle = 45, hjust = 1),                 # angled x labels
    panel.grid = element_blank(),                                     # remove gridlines
    legend.position = "none"                                          # optional: no legend
  )

  xlab("Timepoint") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # centered title
    axis.text.x = element_text(angle = 45, hjust = 1),                 # angled x labels
    panel.grid = element_blank(),                                     # remove gridlines
    legend.position = "none"                                          # optional: no legend
  )





  myeloid_subset <- subset(spatial, idents = "Ctss+ Myeloid Cells")

  cell_counts <- table(myeloid_subset$timepoint)

  time_labels <- paste0(names(cell_counts), "\n", as.vector(cell_counts), " Cells")

  VlnPlot(
    myeloid_subset,
    features = "Gja1",
    group.by = "timepoint",
    pt.size = 0,  # no points
    cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
  ) +
    ggtitle("Gja1 Expression in Ctss+ Myeloid Cells Across Timepoints") +
    ylab("Expression Level") +
    xlab("Timepoint") +
    scale_x_discrete(labels = time_labels) +  # use stacked labels
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # centered title
      axis.text.x = element_text(size = 10, vjust = 1),                  # make space for 2-line labels
      panel.grid = element_blank(),                                     # no gridlines
      legend.position = "none"                                          # hide legend
    )







  cell_types <- c(
    "Activated Satellite Cells",
    "FAPs",
    "Myf6+ Mature Myofiber",
    "Type II Mature Myofibers",
    "Type I Mature Myofibers"
  )

  tp_levels <- if (is.factor(spatial$timepoint)) levels(spatial$timepoint) else unique(spatial$timepoint)

  totals <- data.frame(timepoint = spatial$timepoint) |>
    dplyr::count(timepoint, name = "total_cells") |>
    mutate(timepoint = factor(timepoint, levels = tp_levels))

  sub_obj <- subset(spatial, idents = cell_types)

  counts <- data.frame(
    cell_type = Idents(sub_obj),
    timepoint = sub_obj$timepoint
  ) |>
    dplyr::count(cell_type, timepoint, name = "n") |>
    mutate(timepoint = factor(timepoint, levels = tp_levels)) |>
    tidyr::complete(cell_type, timepoint, fill = list(n = 0))

  counts_frac <- counts |>
    left_join(totals, by = "timepoint") |>
    mutate(percent = ifelse(total_cells > 0, 100 * n / total_cells, NA_real_))

  pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

  ggplot(counts_frac, aes(x = timepoint, y = percent, group = cell_type, color = cell_type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.4) +
    scale_color_manual(values = pal) +
    labs(
      title = "Cell-Type Percentages Across Timepoints",
      x = "Timepoint",
      y = "Percent of Total Cells",
      color = "Cell Type"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      axis.line  = element_line(linewidth = 1.2),
      axis.ticks = element_line(linewidth = 1.0),
      legend.position = "right"
    )
