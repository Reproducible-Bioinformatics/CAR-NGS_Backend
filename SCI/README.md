# sci-RNA-seq Pipeline

## Description

This pipeline is designed for processing sci-RNA-seq data, including FASTQ generation, read alignment, barcode information extraction, UMI counting, and final output generation. The script includes multiple stages of processing, such as poly-A trimming, aligning reads using STAR, filtering and sorting BAM files, counting rRNA reads, generating BED files, and computing various statistics. Additionally, the pipeline provides the capability to generate knee plots and a final UMI count matrix for downstream analysis.

## Pipeline Steps

1. FASTQ Generation: The FASTQ files are extracted from BCL files using bcl2fastq. It also includes read loss calculations.
2. Barcode Information Extraction: Read 1 information (RT well, UMI) is integrated into Read 2, adjusting the read name appropriately.
3. Poly-A Trimming: Poly-A tails are trimmed from the reads using TrimGalore.
4. Read Alignment: Reads are aligned using STAR. The aligned reads are stored in BAM format, and alignment statistics are generated.
rRNA Read Filtering: rRNA reads are filtered, sorted, and counted. BAM files are further processed to remove ambiguously mapped reads.
Gene Assignment: Reads are assigned to genes using BED files, and UMI counts are calculated for each sample.
Duplication Rate Calculation: The duplication rate and the proportion of reads that originate from rRNA are computed.
Knee Plot Generation: Knee plots are generated using R scripts to visualize UMI counts per cell.
UMI Count Matrix: A final UMI count matrix is generated, which includes the total number of UMIs per cell, with cells below the cutoff excluded.

## Output

After running the full pipeline, you will generate several output files:

- `alignment_report_complete.txt`: Alignment statistics.
- `rRNA_report_complete.txt`: rRNA statistics.
- `UMI_report_complete.txt`: UMI counts per sample.
- `final-output/rRNA.and.dup.rate.stats`: Final report of rRNA and duplication rates.
- `final-output/knee-plots/`: Knee plots visualizing UMI distributions.
- `prelim.UMI.count.rollup.gz`: Preliminary UMI count matrix.
- `final-output/cell.annotations`: Cell annotation files.

## function usage: sci_fromfastq()

### Description
`sci_fromfastq` is an R function that converts single-cell indexing (SCI) fastq files into a gene expression matrix. It checks the user's environment (Docker or Sudo) and utilizes Docker containers to execute a data processing pipeline. The function generates the gene expression matrix along with quality control (QC) plots and statistics.

### Parameters
- **group**: A character string indicating the user group. The function supports two options:
  - `"docker"`: Use Docker to run the pipeline.
  - `"sudo"`: Use Sudo for permissions.
  
- **folder**: A character string indicating the path to the working directory. This directory should contain the necessary input files.

- **sample.name**: A character string specifying the name of the experiment or sample.

- **UMI.cutoff**: An integer value representing the minimum number of unique molecular identifiers (UMIs) per cell required to consider the cell valid. Cells with UMI counts below this threshold will be excluded from the final matrix.

### Return Value
The function returns a gene expression matrix along with various quality control plots and statistics.

### Example Usage

```r
# Example usage of sci_fromfastq

sci_fromfastq(
  group = "docker",
  folder = "/20tb/ratto/catcheR/tomatrix/",
  sample.name = "H001AS8",
  UMI.cutoff = 500
)