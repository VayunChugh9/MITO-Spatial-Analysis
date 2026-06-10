library(dplyr)
library(Seurat)
library(tibble)

spatial <- readRDS('data/processed/spatial_with_timepoint.rds')

stopifnot("RNA" %in% Assays(spatial))
DefaultAssay(spatial) <- "RNA"

fap <- subset(spatial, idents = "FAPs")

keep_levels <- c("Control","3d","14d","28d")
stopifnot("timepoint" %in% colnames(fap@meta.data))
fap <- subset(fap, subset = timepoint %in% keep_levels)

fap$timepoint <- factor(fap$timepoint, levels = keep_levels)

fap$active_ident <- Idents(fap)
group_vars <- c("active_ident", "timepoint")

fap_bulk <- AggregateExpression(
  fap,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = group_vars
)

Idents(fap_bulk) <- "timepoint"

de_faps_3d_vs_later <- FindMarkers(
  object  = fap_bulk,
  ident.1 = "3d",
  ident.2 = c("Control","14d","28d"),
  test.use = "DESeq2",
  min.pct = 0.1,
  logfc.threshold = 0
) %>% rownames_to_column("gene")

write.csv(de_faps_3d_vs_later, "data/processed/FAPs_3d_vs_Controls14d28d_DESeq2.csv", row.names = FALSE)

head(de_faps_3d_vs_later, 20)
