# Load necessary libraries
library(MAST)
library(edgeR)
library(Matrix)
library(ggplot2)
library(pheatmap)
setwd("/scratch")
# Function to create a unique output directory
create_output_dir <- function(base_name) {
  dir_name <- paste0(base_name, "_feature_selection")
  i <- 1
  final_dir <- dir_name

  while (dir.exists(final_dir)) {
    final_dir <- paste0(dir_name, "_", i)
    i <- i + 1
  }

  dir.create(final_dir, showWarnings = FALSE)
  return(final_dir)
}

# Function to generate a volcano plot
generate_volcano_plot <- function(results, output_file, log2FC_threshold, pvalue_threshold) {
  results$Significance <- ifelse(results$logFC >= log2FC_threshold & results$PValue <= pvalue_threshold, "Upregulated",
                                 ifelse(results$logFC <= -log2FC_threshold & results$PValue <= pvalue_threshold, "Downregulated", "Not significant"))

  plot <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = Significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")) +
    xlim(c(-max(abs(results$logFC)), max(abs(results$logFC)))) +  # Symmetric X-axis
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
    theme_minimal() +
    theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))

  # Save the plot using ggsave
  ggsave(output_file, plot = plot, width = 8, height = 6)
}

# Function to generate a heatmap for DE genes (removing column names)
generate_heatmap <- function(matrix, genes, output_file) {
  if (length(genes) > 1) {
    # Remove column names (cell names) from the heatmap
    pheatmap(matrix[genes, , drop = FALSE], cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", show_rownames = TRUE, show_colnames = FALSE, filename = output_file)
  } else {
    cat("Not enough DE genes available for the heatmap.\n")
  }
}

# Function to perform feature selection using edgeR and MAST (ANOVA)
perform_feature_selection <- function(matrix_file, clustering_file, threshold=0.8, log2FC_threshold=1, pvalue_threshold=0.05, separator=NULL, genes_file=NULL, barcodes_file=NULL,heatmap=FALSE) {
  # Print input parameters for checking
  cat("Input parameters check:\n")
  cat("Matrix file:", matrix_file, "\n")
  cat("Clustering file:", clustering_file, "\n")
  cat("Stability threshold:", threshold, "\n")
  cat("Log2FC threshold:", log2FC_threshold, "\n")
  cat("P-value threshold:", pvalue_threshold, "\n")

  # Determine if the matrix is sparse or dense based on the file extension
  is_sparse <- grepl("\\.mtx(\\.gz)?$", matrix_file)

  # Load the matrix and clustering data
  if (is_sparse) {
    cat("Matrix type: Sparse\n")
    if (is.null(genes_file) || is.null(barcodes_file)) {
      stop("For a sparse matrix, genes_file and barcodes_file must be provided.")
    }
    if (!file.exists(genes_file) || !file.exists(barcodes_file)) {
      stop("Genes file or Barcodes file missing for sparse matrix.")
    }

    # Load the matrix, genes, and barcodes
    matrix <- as.matrix(readMM(gzfile(matrix_file)))
    gene_names <- read.delim(genes_file, header=FALSE)$V1
    cell_barcodes <- read.delim(barcodes_file, header=FALSE)$V1

    # Assign rownames and colnames to the matrix
    if (length(gene_names) != nrow(matrix)) {
      stop("Number of genes does not match the number of rows in the matrix.")
    }
    if (length(cell_barcodes) != ncol(matrix)) {
      stop("Number of barcodes does not match the number of columns in the matrix.")
    }

    rownames(matrix) <- gene_names
    colnames(matrix) <- cell_barcodes

  } else {
    cat("Matrix type: Dense\n")
    if (!file.exists(matrix_file)) {
      stop("The specified dense matrix file does not exist.")
    }
    # Load the dense matrix
    matrix <- as.matrix(read.table(matrix_file, header=TRUE, sep=separator, row.names=1, quote=""))
  }

  # Load the clustering information
  clustering_output <- read.csv(clustering_file)

  # Ensure the clustering file has the necessary columns
  if (!all(c("Cell", "Cluster", "Stability") %in% colnames(clustering_output))) {
    stop("Clustering file must contain 'Cell', 'Cluster', and 'Stability' columns.")
  }

  # Filter matrix to match clustering file cells
  matrix <- matrix[, clustering_output$Cell]
    print(dim(matrix))

  # Directory setup for output
  base_name <- tools::file_path_sans_ext(basename(matrix_file))
  output_dir <- create_output_dir(base_name)
  anova_dir <- file.path(output_dir, "ANOVA_MAST")
  anovalike_dir <- file.path(output_dir, "AnovaLike")
  pairwise_dir <- file.path(output_dir, "PairwiseComparisons")

  dir.create(anova_dir)
  dir.create(anovalike_dir)
  dir.create(pairwise_dir)

  # Perform ANOVA using MAST
  perform_anova_mast(matrix, clustering_output, anova_dir, log2FC_threshold, pvalue_threshold)

  # Perform feature selection (ANOVA-like) for all cells (tutti contro tutti)
  perform_anova_feature_selection(matrix, clustering_output, anovalike_dir, log2FC_threshold, pvalue_threshold)

  # Perform pairwise cluster comparisons
  perform_pairwise_comparisons(matrix, clustering_output, pairwise_dir, log2FC_threshold, pvalue_threshold)

  # Filter based on stability and perform analyses only for stable cells
  stable_cells <- clustering_output$Cell[clustering_output$Stability >= threshold]
  stable_clusters <- clustering_output$Cluster[clustering_output$Stability >= threshold]

  # Ensure non-empty clusters after filtering
  if (length(stable_cells) > 0 && length(unique(stable_clusters)) > 1) {
    stable_matrix <- matrix[, stable_cells]

    # Directories for stable results
    stable_anova_dir <- file.path(output_dir, "ANOVA_MAST_Stable")
    stable_anovalike_dir <- file.path(output_dir, "AnovaLike_Stable")
    stable_pairwise_dir <- file.path(output_dir, "PairwiseComparisons_Stable")

    dir.create(stable_anova_dir)
    dir.create(stable_anovalike_dir)
    dir.create(stable_pairwise_dir)

    # Perform ANOVA and single comparisons on stable cells
    perform_anova_mast(stable_matrix, clustering_output[clustering_output$Cell %in% stable_cells,], stable_anova_dir, log2FC_threshold, pvalue_threshold)
    perform_anova_feature_selection(stable_matrix, clustering_output[clustering_output$Cell %in% stable_cells,], stable_anovalike_dir, log2FC_threshold, pvalue_threshold)
    perform_pairwise_comparisons(stable_matrix, clustering_output[clustering_output$Cell %in% stable_cells,], stable_pairwise_dir, log2FC_threshold, pvalue_threshold)
  } else {
    cat("No valid clusters after stability filtering.\n")
  }

  cat("Feature selection completed.\n")
}

# Function to perform ANOVA using MAST
perform_anova_mast <- function(matrix, clustering_output, output_dir, log2FC_threshold, pvalue_threshold) {
  # Prepare the SingleCellAssay object for MAST
# Verifica se ci sono duplicati o nomi mancanti nelle colonne o nelle righe
if (any(duplicated(rownames(matrix)))) {
  stop("Duplicate gene names found in the matrix.")
}
if (any(duplicated(colnames(matrix)))) {
  stop("Duplicate cell barcodes found in the matrix.")
}
if (any(is.na(rownames(matrix)))) {
  stop("Missing gene names in the matrix.")
}
if (any(is.na(colnames(matrix)))) {
  stop("Missing cell barcodes in the matrix.")
}

# Controlla se le dimensioni della matrice e il clustering sono coerenti
if (ncol(matrix) != length(clustering_output$Cell)) {
  stop("Number of cells in the matrix and clustering output do not match.")
}

# Debug: stampa alcune righe e colonne della matrice
print(rownames(matrix)[1:5])
print(colnames(matrix)[1:5])

# Prepara il SingleCellAssay per MAST
sca <- FromMatrix(exprsArray = log2(matrix + 1), cData = clustering_output)
  # Add cluster information to the assay
  sca$Cluster <- factor(clustering_output$Cluster)

  # Fit the MAST model
  fit <- zlm(~Cluster, sca)

  # Perform likelihood ratio test (LRT)
  lrt <- lrTest(fit, "Cluster")

  p_values <- lrt[, , "Pr(>Chisq)"]
  genes <- rownames(p_values)

  # Create a data frame with gene names and p-values
  mast_results <- data.frame(Gene = genes, PValue = p_values[, "cont"])

  # Perform log2 fold change calculation
  mast_results$logFC <- log2(rowMeans(as.matrix(matrix) + 1e-5))

  # Save ANOVA results
  write.csv(mast_results, file.path(output_dir, "mast_anova_results.csv"),row.names=FALSE)

  # Generate Volcano plot
  volcano_output_file <- file.path(output_dir, "mast_anova_volcano_plot.png")
  generate_volcano_plot(mast_results, volcano_output_file, log2FC_threshold, pvalue_threshold)

  # Filter DE genes based on p-value and log fold change thresholds for heatmap
  DE_genes <- mast_results$Gene[abs(mast_results$logFC) >= log2FC_threshold & mast_results$PValue <= pvalue_threshold]

  if (length(DE_genes) > 0) {
    # Save DE genes heatmap
    heatmap_output_file <- file.path(output_dir, "mast_anova_heatmap.png")
    if(heatmap){generate_heatmap(matrix, DE_genes, heatmap_output_file)}

    # Save the filtered matrix used for the heatmap
    write.csv(matrix[DE_genes, , drop = FALSE], file.path(output_dir, "mast_anova_heatmap_matrix.csv"))
  } else {
    cat("No DE genes passed the thresholds for ANOVA using MAST.\n")
  }

  cat("ANOVA using MAST completed and results saved to:", output_dir, "\n")
}

# Function to perform ANOVA-like feature selection with edgeR
perform_anova_feature_selection <- function(matrix, clustering_output, output_dir, log2FC_threshold, pvalue_threshold) {
  unique_clusters <- unique(clustering_output$Cluster)

  # Create a directory for each cluster comparison
  for (cluster in unique_clusters) {
    comparison_dir <- file.path(output_dir, paste0("Cluster_", cluster, "_vs_all"))
    dir.create(comparison_dir)

    design <- model.matrix(~0 + as.factor(clustering_output$Cluster))
    colnames(design) <- levels(as.factor(clustering_output$Cluster))

    dge <- DGEList(counts=matrix)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit)

    # Extract top genes per cluster based on ANOVA-like analysis
    anova_results <- topTags(lrt, n=nrow(matrix))$table

    # Save DE results matrix
    write.csv(anova_results, file.path(comparison_dir, "anova_DE_results.csv"))

    # Generate Volcano plot
    volcano_output_file <- file.path(comparison_dir, "anova_volcano_plot.png")
    generate_volcano_plot(anova_results, volcano_output_file, log2FC_threshold, pvalue_threshold)
    # Filter DE genes based on thresholds for heatmap
    DE_genes <- rownames(anova_results[abs(anova_results$logFC) >= log2FC_threshold & anova_results$PValue <= pvalue_threshold,])

    if (length(DE_genes) > 0) {
      # Save DE genes heatmap
      heatmap_output_file <- file.path(comparison_dir, "anova_heatmap.png")

      if(heatmap){generate_heatmap(matrix, DE_genes, heatmap_output_file)}

      # Save the filtered matrix used for the heatmap
      write.csv(matrix[DE_genes, , drop = FALSE], file.path(comparison_dir, "anova_heatmap_matrix.csv"))
    } else {
      cat("No DE genes passed the thresholds for cluster", cluster, "in the ANOVA-like comparison.\n")
    }

    cat("ANOVA-like feature selection for cluster", cluster, "vs all others saved to:", comparison_dir, "\n")
  }
}

# Function to perform pairwise cluster comparisons (e.g., 1 vs 2, 1 vs 3, etc.)
perform_pairwise_comparisons <- function(matrix, clustering_output, output_dir, log2FC_threshold, pvalue_threshold) {
  unique_clusters <- unique(clustering_output$Cluster)

  # Perform all pairwise comparisons
  for (i in 1:(length(unique_clusters))) {
    for (j in 1:length(unique_clusters)) {
    if(i!=j){
      cluster_i <- unique_clusters[i]
      cluster_j <- unique_clusters[j]

      cat("Performing comparison for clusters:", cluster_i, "vs", cluster_j, "\n")

      # Create a directory for each pairwise comparison
      comparison_dir <- file.path(output_dir, paste0("Cluster_", cluster_i, "_vs_", cluster_j))
      dir.create(comparison_dir)

      # Create comparison groups (cluster_i vs cluster_j)
      comparison_cells <- clustering_output$Cell[clustering_output$Cluster %in% c(cluster_i, cluster_j)]
      comparison_matrix <- matrix[, comparison_cells]

      group <- factor(ifelse(clustering_output$Cluster[clustering_output$Cell %in% comparison_cells] == cluster_i, cluster_i, cluster_j))
      design <- model.matrix(~0 + group)
      colnames(design) <- c(cluster_i, cluster_j)

      dge <- DGEList(counts = comparison_matrix)
      dge <- calcNormFactors(dge)
      dge <- estimateDisp(dge, design)
      fit <- glmFit(dge, design)
      lrt <- glmLRT(fit, contrast = c(1, -1))  # Compare cluster_i vs cluster_j

      # Extract top genes for this comparison
      pairwise_results <- topTags(lrt, n=nrow(matrix))$table

      # Save DE results matrix
      write.csv(pairwise_results, file.path(comparison_dir, paste0("DE_results_cluster_", cluster_i, "_vs_", cluster_j, ".csv")))

      # Generate Volcano plot
      volcano_output_file <- file.path(comparison_dir, paste0("volcano_cluster_", cluster_i, "_vs_", cluster_j, ".png"))
      generate_volcano_plot(pairwise_results, volcano_output_file, log2FC_threshold, pvalue_threshold)

      # Filter DE genes based on thresholds for heatmap
      DE_genes <- rownames(pairwise_results[abs(pairwise_results$logFC) >= log2FC_threshold & pairwise_results$PValue <= pvalue_threshold,])

      if (length(DE_genes) > 0) {
        # Save DE genes heatmap
        heatmap_output_file <- file.path(comparison_dir, paste0("heatmap_cluster_", cluster_i, "_vs_", cluster_j, ".png"))
        if(heatmap){generate_heatmap(comparison_matrix, DE_genes, heatmap_output_file)}

        # Save the filtered matrix used for the heatmap
        write.csv(comparison_matrix[DE_genes, , drop = FALSE], file.path(comparison_dir, paste0("heatmap_matrix_cluster_", cluster_i, "_vs_", cluster_j, ".csv")))
      } else {
        cat("No DE genes passed the thresholds for cluster", cluster_i, "vs", cluster_j, "\n")
      }

      cat("Pairwise comparison", cluster_i, "vs", cluster_j, "saved to:", comparison_dir, "\n")

    }
    }
  }
}

# Main execution block
args <- commandArgs(trailingOnly = TRUE)

# Parse input arguments
matrix_file <- args[1]  # Path to the matrix file (sparse matrix or dense matrix file)
clustering_file <- args[2]  # Path to the clustering.output.csv file
threshold <- as.numeric(args[3])  # Stability threshold (e.g., 0.8)
log2FC_threshold <- as.numeric(args[4])  # Log2 Fold Change threshold (e.g., 1)
pvalue_threshold <- as.numeric(args[5])  # P-value threshold (e.g., 0.05)

# Optional arguments for sparse matrix handling
separator <- if (length(args) >= 6) args[6] else NULL
genes_file <- if (length(args) >= 7) args[7] else NULL
barcodes_file <- if (length(args) >= 8) args[8] else NULL
heatmap <- if (length(args) >= 9) as.logical(args[9]) else FALSE  # Convert to logical (default: FALSE)

# Call the feature selection function with the parsed arguments
perform_feature_selection(matrix_file, clustering_file, threshold=threshold, log2FC_threshold=log2FC_threshold, pvalue_threshold=pvalue_threshold, separator=separator, genes_file=genes_file, barcodes_file=barcodes_file,heatmap=heatmap)

