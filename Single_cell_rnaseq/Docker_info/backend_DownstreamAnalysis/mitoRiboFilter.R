# Load necessary libraries
library(Matrix)
library(ggplot2)
setwd("/scratch")
# Function to create a unique output directory
create_output_dir <- function(base_name) {
  dir_name <- paste0(base_name, "_mitoribo")
  i <- 1
  final_dir <- dir_name

  while (dir.exists(final_dir)) {
    final_dir <- paste0(dir_name, "_", i)
    i <- i + 1
  }

  dir.create(final_dir, showWarnings = FALSE)
  return(final_dir)
}

# Function to filter and plot matrix based on ribosomal and mitochondrial gene content
filter_and_plot_matrix <- function(matrix_file, mito_range=c(0, 100), ribo_range=c(0, 100), separator=NULL, genes_file=NULL, barcodes_file=NULL) {

  # Print the passed variables for checking
  cat("Input parameters check:\n")
  cat("Matrix file:", matrix_file, "\n")
  cat("Mitochondrial range:", mito_range, "\n")
  cat("Ribosomal range:", ribo_range, "\n")

  # Determine if the matrix is sparse or dense based on the file extension
  is_sparse <- grepl("\\.mtx(\\.gz)?$", matrix_file)

  if (is_sparse) {
    if (!is.null(separator)) {
      warning("Separator is ignored for sparse matrices.")
    }
    cat("Matrix type: Sparse\n")

    # Automatically find genes_file and barcodes_file if not provided
    if (is.null(genes_file) || is.null(barcodes_file)) {
      dir_path <- dirname(matrix_file)
      if (is.null(genes_file)) {
        genes_file <- Sys.glob(file.path(dir_path, "*features.tsv*"))[1]
      }
      if (is.null(barcodes_file)) {
        barcodes_file <- Sys.glob(file.path(dir_path, "*barcodes.tsv*"))[1]
      }
      if (is.null(genes_file) || is.null(barcodes_file)) {
        stop("For a sparse matrix, genes_file and barcodes_file must be provided or found automatically.")
      }
    }

    cat("Genes file:", genes_file, "\n")
    cat("Barcodes file:", barcodes_file, "\n")

  } else {
    if (is.null(separator)) {
      stop("Separator must be provided for a dense matrix.")
    }
    cat("Matrix type: Dense\n")
    cat("Separator:", separator, "\n")
  }

  if (is_sparse) {
    # Sparse matrix: load the matrix, genes, and barcodes files
    if (!file.exists(matrix_file) || !file.exists(genes_file) || !file.exists(barcodes_file)) {
      stop("One or more required files for the sparse matrix are missing.")
    }

    if (grepl("\\.gz$", matrix_file)) {
      matrix <- readMM(gzfile(matrix_file))
    } else {
      matrix <- readMM(matrix_file)
    }

    # Load gene names assuming they are in the first column (not the second column as originally assumed)
    if (grepl("\\.gz$", genes_file)) {
      gene_names <- read.delim(gzfile(genes_file), header=FALSE)$V1
    } else {
      gene_names <- read.delim(genes_file, header=FALSE)$V1
    }

    if (length(gene_names) != nrow(matrix)) {
      stop("Number of genes does not match the number of rows in the matrix.")
    }

    rownames(matrix) <- gene_names

    # Load barcodes
    if (grepl("\\.gz$", barcodes_file)) {
      cell_barcodes <- read.delim(gzfile(barcodes_file), header=FALSE)$V1
    } else {
      cell_barcodes <- read.delim(barcodes_file, header=FALSE)$V1
    }

    colnames(matrix) <- cell_barcodes

  } else {
    # Dense matrix: load the matrix file directly
    if (!file.exists(matrix_file)) {
      stop("The specified dense matrix file does not exist.")
    }

    matrix <- read.table(matrix_file, header=TRUE, sep=separator,row.names=1,quote="")
    gene_names <- rownames(matrix)
  }

  # Extract the base name from the matrix file
  base_name <- tools::file_path_sans_ext(basename(matrix_file))

  # Create the output directory
  output_dir <- create_output_dir(base_name)

  # Calculate ribosomal and mitochondrial percentages
  rps <- grep("(^|[^a-zA-Z])RPS", toupper(gene_names))
  rpl <- grep("(^|[^a-zA-Z])RPL", toupper(gene_names))
  ribosomal_genes <- unique(c(rps, rpl))
  print(ribosomal_genes)
  mitochondrial_genes <- grep("MT:", toupper(gene_names))
  print(mitochondrial_genes)

  expressed_genes <- colSums(matrix >= 3)
  ribosomal_counts <- colSums(matrix[ribosomal_genes, , drop=FALSE] >= 3)
  mitochondrial_counts <- colSums(matrix[mitochondrial_genes, , drop=FALSE] >= 3)

  ribosomal_percentage <- ifelse(expressed_genes > 0, ribosomal_counts / expressed_genes * 100, 0)
  mitochondrial_percentage <- ifelse(expressed_genes > 0, mitochondrial_counts / expressed_genes * 100, 0)

  # Filter cells based on the specified percentages
  valid_cells <- which(ribosomal_percentage >= ribo_range[1] & ribosomal_percentage <= ribo_range[2] &
                         mitochondrial_percentage >= mito_range[1] & mitochondrial_percentage <= mito_range[2])

  filtered_matrix <- matrix[, valid_cells]

  # Create and save the plot
  plot_data <- data.frame(
    Ribosomal = ribosomal_percentage[valid_cells],
    Mitochondrial = mitochondrial_percentage[valid_cells],
    TotalExpressed = expressed_genes[valid_cells]
  )

  plot <- ggplot(plot_data, aes(x = Ribosomal, y = Mitochondrial, color = TotalExpressed)) +
    geom_point(alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = base_name,  # Use the base name as the plot title
      x = "Percentage of Ribosomal Genes Expressed",
      y = "Percentage of Mitochondrial Genes Expressed",
      color = "Total Expressed Genes"
    ) +
    theme_minimal(base_family = "sans") +  # Use the default sans font
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA)    # White plot background
    ) +
    xlim(0, 100) +  # Limit x-axis to 0-100
    ylim(0, 100)    # Limit y-axis to 0-100

  # Save the plot
  plot_path <- file.path(output_dir, "ribosomal_vs_mitochondrial_plot.png")
  ggsave(plot_path, plot)

  # Save the filtered matrix
  if (is_sparse) {
    # Save the sparse matrix in .mtx, .tsv.gz format
    writeMM(filtered_matrix, file.path(output_dir, "filtered_matrix.mtx"))
    write.table(gene_names, file.path(output_dir, "filtered_features.tsv"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(colnames(filtered_matrix), file.path(output_dir, "filtered_barcodes.tsv"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  } else {
    # Save the dense matrix in .txt or .csv format
    write.table(filtered_matrix, file.path(output_dir, "filtered_matrix.txt"), sep=separator, quote=FALSE, col.names=NA)
  }

  cat("Process completed. Plot and filtered matrix saved in:", output_dir, "\n")
}

# Main execution block
args <- commandArgs(trailingOnly = TRUE)

# Parse input arguments
matrix_file <- args[1]  # Path to the matrix file (sparse matrix or dense matrix file)
mito_min <- as.numeric(args[2])  # Min mitochondrial percentage
mito_max <- as.numeric(args[3])  # Max mitochondrial percentage
ribo_min <- as.numeric(args[4])  # Min ribosomal percentage
ribo_max <- as.numeric(args[5])  # Max ribosomal percentage

# Optional arguments for separator and sparse matrices
separator <- if (length(args) >= 6) args[6] else NULL
genes_file <- if (length(args) >= 7) args[7] else NULL
barcodes_file <- if (length(args) >= 8) args[8] else NULL

# Call the function with the parsed arguments
filter_and_plot_matrix(matrix_file, mito_range=c(mito_min, mito_max), ribo_range=c(ribo_min, ribo_max),
                       separator=separator, genes_file=genes_file, barcodes_file=barcodes_file)
