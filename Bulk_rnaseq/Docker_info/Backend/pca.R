#!/usr/bin/env Rscript

library(ggplot2)
setwd("/scratch")
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript pca_script.R <count_matrix_file> <metadata_file> [count_sep] [meta_sep]")
}

count_matrix_file <- args[1]
metadata_file <- args[2]
count_sep <- ifelse(length(args) >= 3, ifelse(args[3] == "tab", "\t", args[3]), ",")
meta_sep <- ifelse(length(args) >= 4, ifelse(args[4] == "tab", "\t", args[4]), ",")

count_matrix <- read.table(count_matrix_file, sep = count_sep, row.names = 1, header = TRUE)
metadata <- read.table(metadata_file, sep = meta_sep, row.names = 1, header = TRUE)

pca_result <- prcomp(t(log2(count_matrix + 1)))
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)
metadata$Sample <- rownames(metadata)
pca_data <- merge(pca_data, metadata, by = "Sample")

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of RNA-Seq Data",
       x = paste0("PC1: ", round(100 * summary(pca_result)$importance[2, 1], 2), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca_result)$importance[2, 2], 2), "% variance")) +
  theme_minimal() +
  theme(legend.title = element_blank())

output_dir <- "/scratch/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Ottieni il nome della matrice di conteggio senza estensione
file_base_name <- tools::file_path_sans_ext(basename(count_matrix_file))

# Crea il nome del file di output per il plot
output_file <- file.path(output_dir, paste0(file_base_name, "_pca_plot.png"))

# Salva il plot con il nome del file basato sulla matrice di conteggio
ggsave(filename = output_file, plot = p, width = 8, height = 6)
