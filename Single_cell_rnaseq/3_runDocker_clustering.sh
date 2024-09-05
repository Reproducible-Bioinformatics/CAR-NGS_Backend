 
#!/bin/bash

# Get the absolute path of the current directory + ./Data
DATA_DIR=$(pwd)/Data

# Define the dense matrix file and clustering parameters
MATRIX_FILE="DENSE_Filtered.csv"  # Name of the dense matrix file
BOOTSTRAP_PERCENTAGE=0.1  # Percentage of cells to remove in each bootstrap iteration
STABILITY_THRESHOLD=0.8  # Stability threshold (e.g., 0.8)
PERMUTATIONS=10  # Number of permutations
SEPARATOR=","  # Separator for the dense matrix
RESOLUTION=0.8  # Seurat clustering resolution

# Run the Docker container
docker run --rm \
  -v ${DATA_DIR}:/scratch \
  repbioinfo/singlecelldownstream \
  Rscript /home/clustering.R /scratch/${MATRIX_FILE} ${BOOTSTRAP_PERCENTAGE} ${STABILITY_THRESHOLD} ${PERMUTATIONS} ${SEPARATOR} NULL NULL ${RESOLUTION}


# Define the sparse matrix files and clustering parameters
MATRIX_FILE="combined_filtered_matrix_with_sample.mtx"  # Sparse matrix file
GENES_FILE="combined_filtered_with_sample_genes.tsv"  # Gene names file
BARCODES_FILE="combined_filtered_with_sample_barcodes.tsv"  # Barcodes file
BOOTSTRAP_PERCENTAGE=0.1  # Percentage of cells to remove in each bootstrap iteration
STABILITY_THRESHOLD=0.8  # Stability threshold (e.g., 0.8)
PERMUTATIONS=10  # Number of permutations
RESOLUTION=0.8  # Seurat clustering resolution

# Run the Docker container
docker run --rm \
  -v ${DATA_DIR}:/scratch \
  repbioinfo/singlecelldownstream \
  Rscript /home/clustering.R /scratch/${MATRIX_FILE} ${BOOTSTRAP_PERCENTAGE} ${STABILITY_THRESHOLD} ${PERMUTATIONS} NULL /scratch/${GENES_FILE} /scratch/${BARCODES_FILE} ${RESOLUTION}
