spatial <- readRDS("data/processed/spatial_with_timepoint.rds")

cell_barcodes <- colnames(spatial)

timepoint <- ifelse(grepl("_1$", cell_barcodes), "Control",
                    ifelse(grepl("_2$", cell_barcodes), "3d",
                           ifelse(grepl("_3$", cell_barcodes), "14d",
                                  ifelse(grepl("_4$", cell_barcodes), "28d", NA))))

spatial$timepoint <- factor(timepoint, levels = c("Control", "3d", "14d", "28d"))

table(spatial$timepoint, useNA = "ifany")
