library(dplyr)
library(FNN)
library(grid)
library(Matrix)
library(pheatmap)
library(Seurat)
library(spdep)

spatial_3d <- subset(spatial, subset = timepoint == "Control")

meta <- spatial_3d@meta.data
subset_idx <- which(meta$x_centroid >= 5500 & meta$x_centroid <= 10000 &
                      meta$y_centroid >= 0 & meta$y_centroid <= 5000)
spatial_3d <- subset(spatial_3d, cells = rownames(meta)[subset_idx])

gene <- "Gja1"                      # adjust if feature name differs
gexp <- FetchData(spatial_3d, vars = gene)[, 1]

spatial_3d$Gja1_expr <- gexp
spatial_3d$Gja1_pos  <- gexp > 0.25

new_ids <- as.character(Idents(spatial_3d))
new_ids[new_ids == "FAPs" & spatial_3d$Gja1_pos] <- "Gja1+ FAPs"
Idents(spatial_3d) <- factor(new_ids)

if (exists("celltype_colors") && is.character(celltype_colors)) {
  if ("FAPs" %in% names(celltype_colors)) {
    celltype_colors["Gja1+ FAPs"] <- celltype_colors["FAPs"]
  }
}

coords <- as.matrix(spatial_3d@meta.data[, c("x_centroid", "y_centroid")])
cell_types <- Idents(spatial_3d)
types <- levels(cell_types)

k <- 10
nn <- get.knn(coords, k = k)

co_matrix <- matrix(0L, nrow = length(types), ncol = length(types),
                    dimnames = list(types, types))

for (i in seq_len(nrow(coords))) {
  A <- as.character(cell_types[i])
  Btypes <- as.character(cell_types[nn$nn.index[i, ]])
  for (B in Btypes) co_matrix[A, B] <- co_matrix[A, B] + 1L
}

row_sums <- rowSums(co_matrix)
row_sums[row_sums == 0] <- 1
norm_mat <- sweep(co_matrix, 1, row_sums, "/")

pheatmap(
  norm_mat,
  main = "Control ROI #3 Cell-type Co-occurrence",
  color = colorRampPalette(c("white", "steelblue"))(100),
  border_color = NA
)

set.seed(123)

labels <- Idents(spatial_3d)
types  <- levels(labels)
Tn     <- length(types)
A_idx  <- as.integer(labels)
B_idx_mat <- matrix(as.integer(labels)[nn$nn.index], ncol = k)

count_AB <- function(A_idx, B_idx_mat, Tn) {
  out <- matrix(0L, nrow = Tn, ncol = Tn, dimnames = list(types, types))
  for (a in seq_len(Tn)) {
    rows <- which(A_idx == a)
    if (length(rows)) {
      Bs <- as.vector(B_idx_mat[rows, , drop = FALSE])
      out[a, ] <- tabulate(Bs, nbins = Tn)
    }
  }
  out
}

obs_counts <- count_AB(A_idx, B_idx_mat, Tn)
row_totals <- as.numeric(table(factor(A_idx, levels = 1:Tn))) * k
obs_frac   <- sweep(obs_counts, 1, pmax(row_totals, 1), "/")

B <- 1000
perm_ge <- matrix(0L, Tn, Tn)
perm_le <- matrix(0L, Tn, Tn)
perm_mean <- matrix(0, Tn, Tn)

for (b in seq_len(B)) {
  A_perm <- sample(A_idx)
  counts_b <- count_AB(A_perm, B_idx_mat, Tn)
  frac_b <- sweep(counts_b, 1, pmax(row_totals, 1), "/")
  perm_ge <- perm_ge + (frac_b >= obs_frac)
  perm_le <- perm_le + (frac_b <= obs_frac)
  perm_mean <- perm_mean + frac_b
}
perm_mean <- perm_mean / B

p_hi <- (perm_ge + 1) / (B + 1)
p_lo <- (perm_le + 1) / (B + 1)
p_two <- pmin(1, 2 * pmin(p_hi, p_lo))
p_adj <- matrix(p.adjust(as.vector(p_two), method = "BH"),
                nrow = Tn, ncol = Tn, dimnames = dimnames(obs_frac))

alpha <- 0.05
class_code <- matrix(0L, Tn, Tn, dimnames = dimnames(obs_frac))
class_code[p_adj < alpha & obs_frac > perm_mean] <-  1L
class_code[p_adj < alpha & obs_frac < perm_mean] <- -1L

row_clust <- hclust(dist(norm_mat))
col_clust <- hclust(dist(t(norm_mat)))
ord_rows <- rownames(norm_mat)[row_clust$order]
ord_cols <- colnames(norm_mat)[col_clust$order]

class_ord <- class_code[ord_rows, ord_cols, drop = FALSE]

diag(class_ord) <- NA

cols <- c(
  "negative" = "#0072B2",  # blue
  "random"   = "#D7D7D7",  # light neutral gray
  "positive" = "#D55E00"   # vermillion
)
breaks <- c(-1.5, -0.5, 0.5, 1.5)

pheatmap(
  class_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = unname(cols),
  breaks = breaks,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("negative", "random", "positive"),
  na_col = "#F5F5F5",
  border_color = "white",          # subtle gridlines
  cellheight = 14,
  cellwidth = 14,
  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = "45",
  main = "Control ROI #3 Co-occurrence",
  treeheight_row = 0,
  treeheight_col = 0,
  legend = TRUE,
  name = " "                       # visually blank legend title
)



expr <- FetchData(spatial_3d, vars = "Gja1")$Gja1

k <- 10
knn_list <- knearneigh(coords, k = k)
nb <- knn2nb(knn_list)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_result <- moran.test(expr, lw, zero.policy = TRUE)
moran_result
