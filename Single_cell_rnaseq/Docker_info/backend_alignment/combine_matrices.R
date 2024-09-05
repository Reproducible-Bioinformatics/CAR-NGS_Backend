#!/usr/bin/env Rscript

# Carica le librerie necessarie
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

# Funzione per leggere una matrice di conteggio filtrata da una directory
read_matrix <- function(matrix_dir) {
  barcode_path <- file.path(matrix_dir, "barcodes_modified.tsv.gz")
  features_path <- file.path(matrix_dir, "features.tsv.gz")
  matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")

  # Usa on.exit per chiudere automaticamente le connessioni
  barcodes <- readLines(gzfile(barcode_path))
  on.exit(closeAllConnections(), add = TRUE)

  # Leggi il file features.tsv.gz e processa le annotazioni
  features <- read.delim(gzfile(features_path), header = FALSE, stringsAsFactors = FALSE)
  on.exit(closeAllConnections(), add = TRUE)

  # Filtra per "Gene Expression" e combina le colonne 1 e 2
  features <- features[features$V3 == "Gene Expression", ]
  annotated_genes <- paste0(features$V1, ":", features$V2)

  # Leggi la matrice di conteggio
  matrix <- readMM(gzfile(matrix_path))
  on.exit(closeAllConnections(), add = TRUE)

  # Imposta i nomi delle righe come gli ID annotati
  rownames(matrix) <- annotated_genes
  colnames(matrix) <- barcodes

  return(matrix)
}

# Funzione per calcolare le percentuali e generare il plot
calculate_and_plot <- function(matrix, output_dir) {
  gene_names <- rownames(matrix)

  # Identifica i geni ribosomali e mitocondriali
  rps <- grep("RPS", toupper(gene_names))
  rpl <- grep("RPL", toupper(gene_names))
  ribosomal_genes <- unique(c(rps, rpl))
  mitochondrial_genes <- grep("MT:", toupper(gene_names))

  cat("Numero totale di geni ribosomali identificati:", length(ribosomal_genes), "\n")
  cat("Numero totale di geni mitocondriali identificati:", length(mitochondrial_genes), "\n")

  # Calcola il numero di geni espressi per cellula (almeno 3 conte)
  expressed_genes <- colSums(matrix >= 3)


  # Calcola la percentuale di geni ribosomali espressi per cellula
  ribosomal_counts <- colSums(matrix[ribosomal_genes, , drop=FALSE] >= 3)

  ribosomal_percentage <- ifelse(expressed_genes > 0, ribosomal_counts / expressed_genes * 100, 0)

  # Calcola la percentuale di geni mitocondriali espressi per cellula
  mitochondrial_counts <- colSums(matrix[mitochondrial_genes, , drop=FALSE] >= 3)

  mitochondrial_percentage <- ifelse(expressed_genes > 0, mitochondrial_counts / expressed_genes * 100, 0)

  # Numero totale di geni espressi per cellula (almeno 3 conte)
  total_expressed_genes <- expressed_genes
  cat("Distribuzione del numero totale di geni espressi per cellula:\n")

  # Controlla le dimensioni dei vettori
  cat("Dimensioni dei vettori:\n")
  cat("Ribosomal Percentage:", length(ribosomal_percentage), "\n")
  cat("Mitochondrial Percentage:", length(mitochondrial_percentage), "\n")
  cat("Total Expressed Genes:", length(total_expressed_genes), "\n")

  # Creazione del dataframe per il plot
  plot_data <- data.frame(
    Ribosomal = ribosomal_percentage,
    Mitochondrial = mitochondrial_percentage,
    TotalExpressed = total_expressed_genes
  )

  # Genera il plot
  plot <- ggplot(plot_data, aes(x = Ribosomal, y = Mitochondrial, color = TotalExpressed)) +
    geom_point(alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = "Ribosomal vs Mitochondrial Gene Expression",
      x = "Percentage of Ribosomal Genes Expressed",
      y = "Percentage of Mitochondrial Genes Expressed",
      color = "Total Expressed Genes"
    ) +
    theme_minimal()

  # Salva il plot
  plot_path <- file.path(output_dir, "ribosomal_vs_mitochondrial_plot.png")
  ggsave(plot_path, plot)

  cat("Plot salvato in:", plot_path, "\n")
}

# Funzione per salvare barcodes e geni annotati
save_barcodes_and_genes <- function(matrix, dir_path, prefix) {
  # Salva i nomi delle righe (geni annotati)
  writeLines(rownames(matrix), file.path(dir_path, paste0(prefix, "_genes.tsv")))

  # Salva i barcodes
  writeLines(colnames(matrix), file.path(dir_path, paste0(prefix, "_barcodes.tsv")))
}

# Funzione per combinare le matrici di conteggio filtrate da più cartelle
combine_matrices <- function(sample_dirs) {
  all_matrices <- list()
  for (sample_dir in sample_dirs) {
    matrix_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
    matrix <- read_matrix(matrix_dir)

    # Aggiunge la matrice allaF lista
    all_matrices[[sample_dir]] <- matrix
  }

  # Combina tutte le matrici in una sola matrice grande
  combined_matrix <- do.call(cbind, all_matrices)
  return(combined_matrix)
}

# Funzione per salvare una matrice sparsa in formato Matrix Market (.mtx)
save_sparse_matrix <- function(matrix, dir_path, prefix) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  writeMM(matrix, file = file.path(dir_path, paste0(prefix, ".mtx")))
  cat("Matrice salvata in:", file.path(dir_path, paste0(prefix, ".mtx")), "\n")
}

# Elenco degli argomenti passati al comando
args <- commandArgs(trailingOnly = TRUE)

# Tutti gli argomenti sono directory dei campioni
sample_dirs <- args

# Directory base per l'output
output_base_dir <- "/scratch/output"

# Percorso di output per la matrice filtrata combinata
filtered_with_sample_dir <- file.path(output_base_dir, "filtered_with_sample")

# Combina le matrici filtrate
filtered_matrix_with_sample <- combine_matrices(sample_dirs)

# Salva la matrice filtrata combinata in formato Matrix Market (.mtx)
save_sparse_matrix(filtered_matrix_with_sample, filtered_with_sample_dir, "combined_filtered_matrix_with_sample")

# Salva anche i nomi delle cellule e dei geni annotati
save_barcodes_and_genes(filtered_matrix_with_sample, filtered_with_sample_dir, "combined_filtered_with_sample")

# Calcola le percentuali e genera il plot
calculate_and_plot(filtered_matrix_with_sample, filtered_with_sample_dir)

cat("Tutte le matrici sono state salvate e il plot è stato generato con successo nelle directory organizzate sotto /scratch/output.\n")
