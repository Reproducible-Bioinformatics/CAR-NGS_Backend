#!/bin/bash

# Get the absolute path of the current directory + ./Data
DATA_DIR=$(pwd)/Data

# Define the dense matrix file and filtering parameters
MATRIX_FILE="DENSE_Filtered.csv"  # Name of the dense matrix file
MITO_MIN=0  # Minimum mitochondrial percentage
MITO_MAX=10  # Maximum mitochondrial percentage
RIBO_MIN=0  # Minimum ribosomal percentage
RIBO_MAX=20  # Maximum ribosomal percentage
SEPARATOR=","  # Separator for the dense matrix

# Run the Docker container
#FOR DENSE
docker run --rm \
  -v ${DATA_DIR}:/scratch \
  repbioinfo/singlecelldownstream \
  Rscript /home/mitoRiboFilter.R ${MATRIX_FILE} ${MITO_MIN} ${MITO_MAX} ${RIBO_MIN} ${RIBO_MAX} ${SEPARATOR}

MATRIX_FILE="combined_filtered_matrix_with_sample.mtx"  # Sparse matrix file
GENES_FILE="combined_filtered_with_sample_genes.tsv"  # Gene names file
BARCODES_FILE="combined_filtered_with_sample_barcodes.tsv"  # Barcodes file
MITO_MIN=0  # Minimum mitochondrial percentage
MITO_MAX=10  # Maximum mitochondrial percentage
RIBO_MIN=0  # Minimum ribosomal percentage
RIBO_MAX=20  # Maximum ribosomal percentage
#FOR SPARSE
  docker run --rm \
  -v ${DATA_DIR}:/scratch \
  repbioinfo/singlecelldownstream \
  Rscript /home/mitoRiboFilter.R ${MATRIX_FILE} ${MITO_MIN} ${MITO_MAX} ${RIBO_MIN} ${RIBO_MAX} NULL ${GENES_FILE} ${BARCODES_FILE}
