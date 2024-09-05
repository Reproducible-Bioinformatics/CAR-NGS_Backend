 
#!/bin/bash

# Get the absolute paths of the data and genome directories
DATA_DIR=$(pwd)/Data
GENOME_DIR=$(pwd)/Genome

# Set the working directory inside the Docker container to /scratch and /genome for the genome files
WORK_DIR=/scratch
GENOME_WORK_DIR=/genome

# Docker image
DOCKER_IMAGE=repbioinfo/carncellranger2

# Parameters for the script
BAMSAVE="true"  # Whether to save the BAM file (true or false)

# Run the Docker command to execute the alignment and indexing process
docker run --rm -v ${DATA_DIR}:${WORK_DIR} \
  -v ${GENOME_DIR}:${GENOME_WORK_DIR} \
  ${DOCKER_IMAGE} \
  bash /home/index_align.sh ${BAMSAVE}
