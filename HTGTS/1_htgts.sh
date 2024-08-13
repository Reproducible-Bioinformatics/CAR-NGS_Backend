#!/bin/bash

if [ $# -lt 6 ]; then
    echo "Usage: $0 <Data_folder> <fastq1> <fastq2> [<expInfo> <expInfo2>] <outDir>"
    exit 1
fi

data_folder=$1
fastq1=$2
fastq2=$3
expInfo=$4
expInfo2=$5
outDir=$6

if [ -z "$fastq2" ]; then
    # Se è stato fornito solo un fastq, non considerare expInfo2
    docker run -iv "$data_folder":/Data repbioinfo/htgts_pipeline_lts /Algorithm/HTGTS_Full.sh -fastq1 /Data/"$fastq1" -expInfo /Data/"$expInfo" -outDir /Data/"$outDir"
else
    docker run -iv "$data_folder":/Data repbioinfo/htgts_pipeline_lts /Algorithm/HTGTS_Full.sh -fastq1 /Data/"$fastq1" -fastq2 /Data/"$fastq2" -expInfo /Data/"$expInfo" -expInfo2 /Data/"$expInfo2" -outDir /Data/"$outDir"
fi
