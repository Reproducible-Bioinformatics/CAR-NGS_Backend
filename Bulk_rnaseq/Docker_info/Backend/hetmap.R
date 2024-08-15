#!/usr/bin/env Rscript

library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript heatmap_script.R <count_matrix_file> <metadata_file> [count_sep] [meta_sep]")
}
setwd("/scratch")
count_matrix_file <- args[1]
metadata_file <- args[2]
count_sep <- ifelse(length(args) >= 3, ifelse(args[3] == "tab", "\t", args[3]), ",")
meta_sep <- ifelse(length(args) >= 4, ifelse(args[4] == "tab", "\t", args[4]), ",")

count_matrix <- read.table(count_matrix_file, sep = count_sep, row.names = 1, header = TRUE)
metadata <- read.table(metadata_file, sep = meta_sep, row.names = 1, header = TRUE)

output_dir <- "/scratch/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Ottieni il nome base del file della matrice di conteggio
file_base_name <- tools::file_path_sans_ext(basename(count_matrix_file))

# Genera il nome del file per la heatmap
heatmap_file <- file.path(output_dir, paste0(file_base_name, "_heatmap_filtered.png"))
if(length(which(rowSums(count_matrix)==0))!=0){
count_matrix=count_matrix[-which(rowSums(count_matrix)==0),]
}
pheatmap(
  mat = count_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = metadata[, "Group", drop = FALSE],
  main = "Heatmap of Significant Genes",
  filename = heatmap_file
)
