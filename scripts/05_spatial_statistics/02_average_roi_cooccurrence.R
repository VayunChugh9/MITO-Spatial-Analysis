library(dplyr)
library(FNN)
library(grid)
library(Matrix)
library(pheatmap)
library(Seurat)
library(spdep)
library(stats)

row_normalize <- function(m) {
  rs <- rowSums(m); rs[rs == 0] <- 1
  sweep(m, 1, rs, "/")
}

align_by_names <- function(m, target_rownames, target_colnames, fill = 0) {
  out <- matrix(fill, nrow = length(target_rownames), ncol = length(target_colnames),
                dimnames = list(target_rownames, target_colnames))
  rn <- intersect(rownames(m), target_rownames)
  cn <- intersect(colnames(m), target_colnames)
  out[rn, cn] <- m[rn, cn, drop = FALSE]
  out
}

count_AB <- function(A_idx, B_idx_mat, types, Tn) {
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

fisher_combine <- function(mats) {
  mats <- lapply(mats, function(M) pmax(M, .Machine$double.xmin))
  S <- Reduce("+", lapply(mats, function(M) -2 * log(M)))
  pchisq(S, df = 2 * length(mats), lower.tail = FALSE)
}

roi_list <- list(
  "Replicate 1" = list(x = c(2500, 6000), y = c(5000, 8000)),
  "Replicate 2" = list(x = c(0, 2500),    y = c(0, 5000)),
  "Replicate 3" = list(x = c(3000, 8000), y = c(1000, 5000))
)

spatial_base <- subset(spatial, subset = timepoint == "3d")

analyze_roi <- function(obj, roi_name, x_rng, y_rng, gene = "Gja1", k = 10, B = 1000, alpha = 0.05) {
  meta <- obj@meta.data
  idx <- which(meta$x_centroid >= x_rng[1] & meta$x_centroid <= x_rng[2] &
                 meta$y_centroid >= y_rng[1] & meta$y_centroid <= y_rng[2])
  if (!length(idx)) stop(sprintf("ROI '%s' had no cells.", roi_name))
  roi <- subset(obj, cells = rownames(meta)[idx])

  gexp <- FetchData(roi, vars = gene)[,1]
  roi$Gja1_expr <- gexp; roi$Gja1_pos <- gexp > 0.25
  new_ids <- as.character(Idents(roi))
  new_ids[new_ids == "FAPs" & roi$Gja1_pos] <- "Gja1+ FAPs"
  Idents(roi) <- factor(new_ids)

  if (exists("celltype_colors") && is.character(celltype_colors)) {
    if ("FAPs" %in% names(celltype_colors)) celltype_colors["Gja1+ FAPs"] <- celltype_colors["FAPs"]
  }

  coords <- as.matrix(roi@meta.data[, c("x_centroid", "y_centroid")])
  labels <- Idents(roi)
  types <- levels(labels); Tn <- length(types)
  nn <- get.knn(coords, k = k)

  co_matrix <- matrix(0L, Tn, Tn, dimnames = list(types, types))
  for (i in seq_len(nrow(coords))) {
    A <- as.character(labels[i])
    Btypes <- as.character(labels[nn$nn.index[i, ]])
    for (Btype in Btypes) co_matrix[A, Btype] <- co_matrix[A, Btype] + 1L
  }
  norm_mat <- row_normalize(co_matrix)

  A_idx <- as.integer(labels)
  B_idx_mat <- matrix(as.integer(labels)[nn$nn.index], ncol = k)
  obs_counts <- count_AB(A_idx, B_idx_mat, types, Tn)
  row_totals <- as.numeric(table(factor(A_idx, levels = 1:Tn))) * k
  obs_frac   <- sweep(obs_counts, 1, pmax(row_totals, 1), "/")

  set.seed(123)
  perm_ge <- matrix(0L, Tn, Tn); perm_le <- matrix(0L, Tn, Tn); perm_mean <- matrix(0, Tn, Tn)
  for (b in seq_len(B)) {
    A_perm  <- sample(A_idx)
    countsb <- count_AB(A_perm, B_idx_mat, types, Tn)
    fracb   <- sweep(countsb, 1, pmax(row_totals, 1), "/")
    perm_ge   <- perm_ge   + (fracb >= obs_frac)
    perm_le   <- perm_le   + (fracb <= obs_frac)
    perm_mean <- perm_mean + fracb
  }
  perm_mean <- perm_mean / B

  p_hi  <- (perm_ge + 1) / (B + 1)
  p_lo  <- (perm_le + 1) / (B + 1)
  p_two <- pmin(1, 2 * pmin(p_hi, p_lo))
  p_adj <- matrix(p.adjust(as.vector(p_two), method = "BH"),
                  nrow = Tn, ncol = Tn, dimnames = dimnames(obs_frac))

  class_code <- matrix(0L, Tn, Tn, dimnames = dimnames(obs_frac))
  class_code[p_adj < alpha & obs_frac > perm_mean] <-  1L
  class_code[p_adj < alpha & obs_frac < perm_mean] <- -1L
  diag(class_code) <- 0L

  knn_list <- knearneigh(coords, k = k)
  nb <- knn2nb(knn_list)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  moran_res <- suppressWarnings(moran.test(gexp, lw, zero.policy = TRUE))

  list(
    roi_name   = roi_name,
    types      = types,
    norm_mat   = norm_mat,
    obs_frac   = obs_frac,
    perm_mean  = perm_mean,
    class_code = class_code,
    enrich     = obs_frac - perm_mean,
    p_hi       = p_hi,
    p_lo       = p_lo,
    p_two      = p_two,
    p_adj      = p_adj,
    moran      = moran_res
  )
}

roi_results <- lapply(names(roi_list), function(nm) {
  xr <- roi_list[[nm]]$x; yr <- roi_list[[nm]]$y
  analyze_roi(spatial_base, nm, xr, yr, gene = "Gja1", k = 10, B = 1000, alpha = 0.05)
})
names(roi_results) <- names(roi_list)

all_types <- sort(unique(unlist(lapply(roi_results, function(x) x$types))))
get_aligned <- function(res, slot) align_by_names(res[[slot]], all_types, all_types, fill = 0)

norm_list    <- lapply(roi_results, get_aligned, "norm_mat")
enrich_list  <- lapply(roi_results, get_aligned, "enrich")
class_list   <- lapply(roi_results, get_aligned, "class_code")
p_hi_list    <- lapply(roi_results, get_aligned, "p_hi")
p_lo_list    <- lapply(roi_results, get_aligned, "p_lo")

avg_norm   <- Reduce("+", norm_list)   / length(norm_list)
avg_enrich <- Reduce("+", enrich_list) / length(enrich_list)

class_sum <- Reduce("+", class_list)
thresh <- ceiling(length(class_list)/2)
maj_class <- matrix(0L, nrow = length(all_types), ncol = length(all_types),
                    dimnames = list(all_types, all_types))
maj_class[class_sum >=  thresh] <-  1L
maj_class[class_sum <= -thresh] <- -1L

p_comb_pos <- fisher_combine(p_hi_list)
p_comb_neg <- fisher_combine(p_lo_list)
p_comb_two <- pmin(1, 2 * pmin(p_comb_pos, p_comb_neg))
p_comb_adj <- matrix(p.adjust(as.vector(p_comb_two), method = "BH"),
                     nrow = length(all_types), ncol = length(all_types),
                     dimnames = list(all_types, all_types))

dir_mat <- matrix(0L, nrow = length(all_types), ncol = length(all_types),
                  dimnames = list(all_types, all_types))
dir_mat[p_comb_pos < p_comb_neg] <-  1L
dir_mat[p_comb_pos > p_comb_neg] <- -1L

comb_class <- matrix(0L, nrow = length(all_types), ncol = length(all_types),
                     dimnames = list(all_types, all_types))
alpha <- 0.05
comb_class[p_comb_adj < alpha & dir_mat ==  1L] <-  1L
comb_class[p_comb_adj < alpha & dir_mat == -1L] <- -1L

sym_mat <- (avg_norm + t(avg_norm)) / 2
ord <- rownames(sym_mat)[hclust(dist(sym_mat))$order]   # or simply: ord <- all_types

avg_norm_ord   <- avg_norm  [ord, ord, drop = FALSE]
avg_enrich_ord <- avg_enrich[ord, ord, drop = FALSE]
maj_class_ord  <- maj_class [ord, ord, drop = FALSE]
comb_class_ord <- comb_class[ord, ord, drop = FALSE]
p_comb_adj_ord <- p_comb_adj[ord, ord, drop = FALSE]

diag(maj_class_ord)  <- NA_integer_
diag(comb_class_ord) <- NA_integer_

cols_disc <- c("negative" = "#0072B2", "random" = "#D7D7D7", "positive" = "#D55E00")

pheatmap(
  avg_norm_ord,
  main = "3d — Mean co-occurrence (row-normalized) across ROIs",
  color = colorRampPalette(c("white", "steelblue"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA
)

enrich_lim <- max(abs(range(avg_enrich_ord, finite = TRUE)))
pheatmap(
  avg_enrich_ord,
  main = "3d — Mean enrichment score (obs − perm_mean) across ROIs",
  color = colorRampPalette(c("#0072B2", "white", "#D55E00"))(101),
  breaks = seq(-enrich_lim, enrich_lim, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA
)

pheatmap(
  maj_class_ord,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = unname(cols_disc),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("negative", "ns", "positive"),
  na_col = "white",                 # diagonal shows as white
  border_color = "white",
  cellheight = 14, cellwidth = 14, fontsize = 9, angle_col = "45",
  main = "3d — Majority-vote co-occurrence (diagonal masked)"
)

pheatmap(
  comb_class_ord,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = unname(cols_disc),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("negative", "ns", "positive"),
  na_col = "white",                 # diagonal shows as white
  border_color = "white",
  cellheight = 14, cellwidth = 14, fontsize = 9, angle_col = "45",
  main = "3d Co-occurence by Cell Type"
)

neglog10 <- function(x) { y <- -log10(pmax(x, .Machine$double.xmin)); y[!is.finite(y)] <- 0; y }
pheatmap(
  neglog10(p_comb_adj_ord),
  main = "3d Co-occurence by Cell Type)",
  color = colorRampPalette(c("white", "black"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA
)

moran_table <- do.call(rbind, lapply(roi_results, function(x) {
  est <- tryCatch(unname(x$moran$estimate[["Moran I statistic standard deviate"]]), error = function(e) NA_real_)
  if (is.na(est)) est <- tryCatch(unname(x$moran$estimate[[1]]), error = function(e) NA_real_)
  cbind(ROI = x$roi_name, I = est, p = x$moran$p.value)
}))
print(moran_table)
cat(sprintf("\n3d — Mean Moran's I summary across ROIs: %.3f\n",
            mean(suppressWarnings(as.numeric(moran_table[, "I"])), na.rm = TRUE)))
