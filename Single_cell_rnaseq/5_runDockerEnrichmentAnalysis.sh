#!/bin/bash

# Get the absolute path of the current directory (assumed to be where the data is stored)
DATA_DIR=$(pwd)/Data

# Set the working directory inside the Docker container to /scratch
WORK_DIR=/scratch

# Docker image
DOCKER_IMAGE=repbioinfo/singlecelldownstream

# Input parameters (hardcoded in this case)
FILE_NAME="anova_DE_results.csv"  # File name for DE results
SPECIES="dmelanogaster"           # Species
SOURCE="KEGG"                     # Source for enrichment analysis
SEPARATOR=","                     # Separator for the file ("," in this case)
MAX_TERMS=20                      # Maximum number of terms to show in the enrichment plot

# Docker command to run the enrichment analysis
docker run --rm -v ${DATA_DIR}:${WORK_DIR} \
  ${DOCKER_IMAGE} \
  Rscript /home/enrichment_analysis.R ${FILE_NAME} ${SPECIES} ${SOURCE} ${SEPARATOR} ${MAX_TERMS}
