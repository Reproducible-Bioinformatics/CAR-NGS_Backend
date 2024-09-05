#!/bin/bash

# Get the absolute path of the current directory + ./Data
DATA_DIR=$(pwd)/Data

# Define the sparse matrix files and clustering parameters
MATRIX_FILE="combined_filtered_matrix_with_sample.mtx"  # Sparse matrix file
CLUSTERING_FILE="combined_filtered_matrix_with_sample_clustering_stability.output.csv"  # Clustering file
GENES_FILE="combined_filtered_with_sample_genes.tsv"  # Gene names file
BARCODES_FILE="combined_filtered_with_sample_barcodes.tsv"  # Barcodes file
THRESHOLD=0  # Stability threshold
LOG2FC_THRESHOLD=1  # Log2 Fold Change threshold
PVALUE_THRESHOLD=0.05  # P-value threshold
HEATMAP="FALSE"  # Heatmap option

# Run the Docker container
#docker run --rm \
#  -v ${DATA_DIR}:/scratch \
#  repbioinfo/singlecelldownstream \
#  Rscript /home/featureSelection.R ${MATRIX_FILE} ${CLUSTERING_FILE} ${THRESHOLD} ${LOG2FC_THRESHOLD} ${PVALUE_THRESHOLD} NULL ${GENES_FILE} ${BARCODES_FILE} ${HEATMAP}

# Define the dense matrix file and clustering parameters
MATRIX_FILE="DENSE_Filtered.csv"  # Dense matrix file
CLUSTERING_FILE="DENSE_Filtered_clustering_stability.output.csv"  # Clustering file
SEPARATOR=","  # Separator for the dense file
THRESHOLD=0  # Stability threshold
LOG2FC_THRESHOLD=1  # Log2 Fold Change threshold
PVALUE_THRESHOLD=0.05  # P-value threshold
HEATMAP="FALSE"  # Heatmap option

# Run the Docker container
docker run --rm \
  -v ${DATA_DIR}:/scratch \
  repbioinfo/singlecelldownstream \
  Rscript /home/featureSelection.R ${MATRIX_FILE} ${CLUSTERING_FILE} ${THRESHOLD} ${LOG2FC_THRESHOLD} ${PVALUE_THRESHOLD} ${SEPARATOR} NULL NULL ${HEATMAP}
