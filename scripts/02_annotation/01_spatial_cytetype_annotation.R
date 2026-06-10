library(CyteTypeR)
library(dplyr)
library(Seurat)

group_key <- "SCT_snn_res.0.3"
coordinates_key <- "umap"

assay_for_markers <- "SCT"
message("Using assay for markers: ", assay_for_markers)
message("Genes in assay: ", nrow(spatial[[assay_for_markers]]))

Idents(spatial) <- group_key

set.seed(999)
max_n <- 1000

cells_use <- unlist(lapply(levels(Idents(spatial)), function(cl) {
  cells <- WhichCells(spatial, idents = cl)
  if (length(cells) > max_n) sample(cells, max_n) else cells
}), use.names = FALSE)

spatial_ss <- subset(spatial, cells = cells_use)

DefaultAssay(spatial_ss) <- assay_for_markers
Idents(spatial_ss) <- group_key

markers <- FindAllMarkers(
  spatial_ss,
  assay = assay_for_markers,
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.05,
  return.thresh = 1,
  recorrect_umi = FALSE
)

if (!"gene" %in% colnames(markers)) {
  markers$gene <- rownames(markers)
}

if (!"cluster" %in% colnames(markers) && "ident" %in% colnames(markers)) {
  markers <- dplyr::rename(markers, cluster = ident)
}
stopifnot("cluster" %in% colnames(markers))

markers$cluster <- as.character(markers$cluster)


blacklist_patterns <- c(
  "^Tr[abdg][vcj]", "^Igh", "^Igk", "^Igl",
  "^Rpl", "^Rps", "^mt-",
  "^Hba", "^Hbb",
  "^Malat1$", "^Neat1$", "^Xist$"
)
blacklist_regex <- paste(blacklist_patterns, collapse="|")

markers_filt <- markers %>%
  filter(!grepl(blacklist_regex, gene, ignore.case = TRUE))

message("Markers rows after filtering: ", nrow(markers_filt))
message("Clusters represented in markers: ", length(unique(markers_filt$cluster)))

prepped_data <- PrepareCyteTypeR(
  spatial,
  markers_filt,
  n_top_genes = 50,
  group_key = group_key,
  aggregate_metadata = TRUE,
  coordinates_key = coordinates_key
)

study_context <- "
Single-cell transcriptomics data from mouse skeletal muscle tissue collected across a regeneration time course.
The samples include uninjured control muscle and injured muscle harvested at 3, 14, and 28 days post injury.
The study focuses on injury-induced transcriptional programs, muscle development pathways, inflammatory responses, and stromal remodeling during skeletal muscle regeneration.
"

metadata <- list(
  title = "Mouse skeletal muscle regeneration scRNA-seq",
  run_label = "cytetype_seurat_SCT_snn_res.0.3_v1",
  experiment_name = "muscle_regen_timecourse_mouse"
)

spatial_annot <- CyteTypeR(
  obj = spatial,
  prepped_data = prepped_data,
  study_context = study_context,
  metadata = metadata
)

grep("^cytetype", colnames(spatial_annot@meta.data), value = TRUE)

DimPlot(
  spatial_annot,
  reduction = coordinates_key,
  group.by = paste0("cytetype_", group_key),
  label = TRUE
)

table(
  cyte = spatial_annot[[paste0("cytetype_", group_key)]][,1],
  manual = spatial_annot$manual_ident
)
