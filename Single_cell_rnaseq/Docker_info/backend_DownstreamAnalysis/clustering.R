# Load necessary libraries
library(Matrix)
library(ggplot2)
library(future)
library(future.apply)
setwd("/scratch")
# Function to create a unique output directory
create_output_dir <- function(base_name) {
  dir_name <- paste0(base_name, "_clustering")
  i <- 1
  final_dir <- dir_name

  while (dir.exists(final_dir)) {
    final_dir <- paste0(dir_name, "_", i)
    i <- i + 1
  }

  dir.create(final_dir, showWarnings = FALSE)
  return(final_dir)
}

# Function to calculate the Jaccard index between two sets of cells
calculate_jaccard_index <- function(original_cells, new_cells) {
  intersect_size <- length(intersect(original_cells, new_cells))
  union_size <- length(union(original_cells, new_cells))
  return(intersect_size / union_size)
}

# Generalized custom clustering function using Seurat, taking resolution as an explicit parameter
custom_clustering <- function(expression_matrix, resolution) {
  library(Seurat)
  # Example using Seurat (you can replace this with another algorithm)
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  seurat_object <- FindNeighbors(seurat_object)
  seurat_object <- FindClusters(seurat_object, resolution = resolution)
  seurat_object <- RunUMAP(seurat_object, dims = 1:10)  # You can modify the dimensionality reduction here

  # Return a named vector of clusters and UMAP coordinates
  clusters <- seurat_object$seurat_clusters
  umap_coords <- Embeddings(seurat_object, "umap")

  return(list(clusters = clusters, umap = umap_coords))
}

# Function to perform clustering and stability analysis
perform_clustering_and_stability <- function(matrix_file, separator=NULL, genes_file=NULL, barcodes_file=NULL, bootstrap_percentage=0.1, stability_threshold=0.8, permutations=100, resolution=0.8) {

  # Print the passed variables for checking
  cat("Input parameters check:\n")
  cat("Matrix file:", matrix_file, "\n")
  cat("Bootstrap percentage:", bootstrap_percentage, "\n")
  cat("Stability threshold:", stability_threshold, "\n")
  cat("Permutations:", permutations, "\n")
  cat("Resolution:", resolution, "\n")

  # Determine if the matrix is sparse or dense based on the file extension
  is_sparse <- grepl("\\.mtx(\\.gz)?$", matrix_file)

  if (is_sparse) {
    if (!is.null(separator)) {
      warning("Separator is ignored for sparse matrices.")
    }
    cat("Matrix type: Sparse\n")

    # Ensure that genes_file and barcodes_file are provided
    if (is.null(genes_file) || is.null(barcodes_file)) {
      stop("For a sparse matrix, genes_file and barcodes_file must be provided.")
    }

    # Load the sparse matrix
    matrix <- readMM(gzfile(matrix_file))
    gene_names <- read.delim(gzfile(genes_file), header=FALSE)$V1
    cell_barcodes <- read.delim(gzfile(barcodes_file), header=FALSE)$V1
    rownames(matrix) <- gene_names
    colnames(matrix) <- cell_barcodes

  } else {
    # Dense matrix: load the matrix file directly
    if (!file.exists(matrix_file)) {
      stop("The specified dense matrix file does not exist.")
    }

    matrix <- read.table(matrix_file, header=TRUE, sep=separator, row.names=1, quote="")
  }

  # Apply the custom clustering function and check its output
  clustering_result <- custom_clustering(matrix, resolution)
  clusters <- clustering_result$clusters
  umap_coords <- clustering_result$umap

  # Create output directory
  base_name <- tools::file_path_sans_ext(basename(matrix_file))
  output_dir <- create_output_dir(base_name)

  # Save clustering results with stability to CSV
  clustering_output <- data.frame(Cell = names(clusters), Cluster = clusters, UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2])

  # Perform stability analysis
  stability_scores <- numeric(ncol(matrix))
  names(stability_scores) <- colnames(matrix)

  plan("multisession", workers = availableCores())  # Use parallel processing

  # Parallel bootstrap loop
  bootstrap_results <- future_lapply(1:permutations, function(i) {
    # Sample 90% of the cells
    sampled_cells <- sample(colnames(matrix), size = floor((1 - bootstrap_percentage) * ncol(matrix)))
    sampled_matrix <- matrix[, sampled_cells]

    # Apply clustering function to the sampled dataset
    sampled_clusters <- custom_clustering(sampled_matrix, resolution)$clusters

    # Initialize Jaccard score for this bootstrap
    jaccard_scores <- numeric(length(clusters))
    names(jaccard_scores) <- names(clusters)

    # Calculate Jaccard scores for each cluster
    for (cluster in unique(clusters)) {
      original_cells <- names(clusters[clusters == cluster])

      # Identify the matching cluster in the new clustering based on Jaccard index
      best_jaccard <- 0
      best_cluster <- NULL

      for (new_cluster in unique(sampled_clusters)) {
        new_cells <- names(sampled_clusters[sampled_clusters == new_cluster])

        jaccard_index <- calculate_jaccard_index(original_cells, new_cells)

        if (jaccard_index > best_jaccard) {
          best_jaccard <- jaccard_index
          best_cluster <- new_cluster
        }
      }

      # If the best Jaccard index is greater than or equal to the stability threshold, add 1 to all cells in this cluster
      if (best_jaccard >= stability_threshold) {
        matching_cells <- names(sampled_clusters[sampled_clusters == best_cluster])
        jaccard_scores[matching_cells] <- jaccard_scores[matching_cells] + 1
      }
    }

    return(jaccard_scores)
  }, future.seed = TRUE)

  # Combine stability results
  stability_scores <- Reduce("+", bootstrap_results)

  # Calculate final stability percentages
  stability_percentages <- stability_scores / permutations

  # Add stability percentages to clustering output
  clustering_output$Stability <- stability_percentages

  # Save updated clustering output with stability scores
  clustering_output_file_stability <- file.path(output_dir, paste0(base_name, "_clustering_stability.output.csv"))
  write.table(clustering_output, clustering_output_file_stability, sep=",", row.names=FALSE, quote=FALSE)

  cat("Clustering results with stability scores saved to:", clustering_output_file_stability, "\n")

# Save and generate the clustering plots using ggsave
umap_plot_file_clusters <- file.path(output_dir, paste0(base_name, "_umap_clusters.png"))
umap_plot_file_stability <- file.path(output_dir, paste0(base_name, "_umap_stability.png"))

# Plot 1: UMAP colored by cluster
umap_plot_clusters <- ggplot(clustering_output, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(alpha = 0.6) +
  labs(title = "UMAP - Clusters", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # White panel background
        plot.background = element_rect(fill = "white", color = NA))    # White plot background

# Save the first plot
ggsave(umap_plot_file_clusters, plot = umap_plot_clusters, width = 8, height = 6)

# Plot 2: UMAP colored by stability with fixed range 0 to 1
umap_plot_stability <- ggplot(clustering_output, aes(x = UMAP_1, y = UMAP_2, color = Stability)) +
  geom_point(alpha = 0.6) +
  scale_color_gradient(low = "blue", high = "red", limits = c(0, 1)) +  # Fixed range from 0 to 1
  labs(title = "UMAP - Stability", x = "UMAP 1", y = "UMAP 2", color = "Stability") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # White panel background
        plot.background = element_rect(fill = "white", color = NA))    # White plot background

# Save the second plot
ggsave(umap_plot_file_stability, plot = umap_plot_stability, width = 8, height = 6)

cat("Plots saved as PNG files in:", output_dir, "\n")

}

# Main execution block

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)

# Mandatory arguments
matrix_file <- args[1]  # Path to the matrix file (sparse matrix or dense matrix file)
bootstrap_percentage <- as.numeric(args[2])  # Percentage of cells to remove in each bootstrap iteration
stability_threshold <- as.numeric(args[3])  # Stability threshold (e.g., 0.8)
permutations <- as.numeric(args[4])  # Number of permutations

# Optional arguments for data input format
separator <- if (length(args) >= 5) args[5] else NULL  # Separator for dense matrix, ignored for sparse
genes_file <- if (length(args) >= 6) args[6] else NULL  # Genes file (mandatory for sparse matrix)
barcodes_file <- if (length(args) >= 7) args[7] else NULL  # Barcodes file (mandatory for sparse matrix)

# Clustering-related arguments (for SEURAT)
resolution <- if (length(args) >= 8) as.numeric(args[8]) else 0.8  # Resolution parameter for Seurat clustering (default = 0.8)

# If using a different clustering algorithm, replace or append here:
# - For example, if your clustering method requires an additional "k" parameter:
# k <- if (length(args) >= 9) as.numeric(args[9]) else 10  # Example parameter 'k'

# Call the clustering and stability function
perform_clustering_and_stability(matrix_file = matrix_file,
                                 separator = separator,
                                 genes_file = genes_file,
                                 barcodes_file = barcodes_file,
                                 bootstrap_percentage = bootstrap_percentage,
                                 stability_threshold = stability_threshold,
                                 permutations = permutations,
                                 resolution = resolution)

# Note: For other clustering methods:
# If you want to use a different clustering method (e.g., hierarchical clustering, k-means),
# append any additional arguments for that method here and modify the custom_clustering function accordingly.
# Example: If using k-means clustering, you might need to pass 'centers' instead of 'resolution'.
# Example:
# centers <- if (length(args) >= 9) as.numeric(args[9]) else 3  # Number of clusters for k-means

# Add more parameters to the function call above if necessary, depending on the clustering method being used.


