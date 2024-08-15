# RNA-Seq Analysis Pipeline for Drosophila Testes

This repository contains scripts and data to perform a bulk RNA-Seq analysis pipeline for Drosophila testes. The data used in this analysis come from the study "Effect of deletion of CPES on gene expression during spermatogenesis in Drosophila testes" (BioProject: PRJNA860566, GEO: GSE208655). The analysis pipeline includes steps such as indexing and alignment, principal component analysis (PCA), differential expression analysis using DESeq2, and heatmap generation.

## Repository Structure

- `cdsa.sh`: Shell script to execute the complete downstream analysis pipeline.
- `deseq2.sh`: Shell script to execute the DESeq2 analysis.
- `Docker_info`: Directory containing information about the Docker environment used for the analysis.
- `Genome`: Directory containing the genome files split into smaller parts for easy handling. It also contains scripts to recombine these files.
- `heatmap.sh`: Shell script to generate heatmaps from the processed data.
- `index_align.sh`: Shell script to perform genome indexing and alignment of RNA-Seq reads.
- `pca.sh`: Shell script to perform principal component analysis on the RNA-Seq data.
- `README.md`: This file, providing an overview of the repository and instructions for use.

### Data Directory (`Data`)

The `Data` directory contains example datasets to run the pipeline:

- **FASTQ Files**: Raw RNA-Seq data for wild type (wt), CPES mutants (cpes), and rescue samples. These files are compressed and named as `wt1_1.fastq.gz`, `wt1_2.fastq.gz`, `cpes1_1.fastq.gz`, etc.
- **Covariatesstat.csv**: Metadata file containing sample information.
- **gene_count_matrix.csv**: Count matrix used for DESeq2 analysis.
- **filtered_count_matrix.csv**: Filtered count matrix used for downstream analysis.

### Genome Directory (`Genome`)

The `Genome` directory contains the Drosophila genome split into smaller files, which are recombined using the provided scripts:

- `part_aa.gz`, `part_ab.gz`, ...: Parts of a genome file compressed and ready for recombination.
- `split_gtf_01.gz`, `split_gtf_02.gz`, ...: Split parts of the GTF file.
- `extractAndRecombine.sh`: Script to recombine the split GTF and FASTA files into their original format.

## Prerequisites

To run the analysis pipeline, you need the following:

- **Docker**: The analysis pipeline is designed to run inside a Docker container. Ensure Docker is installed and running on your system.

## Prepare the Environment

Ensure Docker is installed and that you have sufficient resources allocated to Docker for the analysis (CPU, memory, etc.). The Docker container used in this repository includes the following dependencies:

- **STAR**: For genome indexing and RNA-Seq read alignment.
- **cutadapt**: For adapter trimming in FASTQ files.
- **DESeq2**: For differential expression analysis.
- **pheatmap**: For generating heatmaps.
- **ggplot2**: For creating PCA plots and other visualizations.
- **VennDiagram**: For generating Venn diagrams of significant genes.
- **AnnotationDbi**: For gene annotation.

All these dependencies are pre-installed in the Docker image `repbioinfo/rnaseqstar_v2`.

## Installation and Setup

1. **Clone the repository**:
   ```bashgit clone https://github.com/yourusername/rnaseq_pipeline.git
   cd rnaseq_pipeline
```
2. **Run the Index and Alignment Script**:
   ```bash./index_align.sh
```
3. **Run the PCA Script**:
   ```bash./pca.sh
```
4. **Run the DESeq2 Script**:
   ```bash./deseq2.sh
```
5. **Generate Heatmaps**:
   ```bash./heatmap.sh
```
6. **Run the Complete Downstream Analysis**:
   ```bash./cdsa.sh
```
## Script Description

### 1. `index_align.sh`

**Description**: This script performs the alignment of RNA-Seq reads to the reference genome using STAR and generates a gene count matrix.

**Docker Input**:
- `/scratch/`: Directory containing FASTQ files (`*.fastq.gz`) of RNA-Seq reads, which are mounted from the `Data` directory on the host.
- `/genome/`: Directory containing the reference genome FASTA file (`*.fa` or `*.fasta`) and GTF annotation file (`*.gtf`), mounted from the `Genome` directory on the host.

**Docker Output**:
- `/scratch/output/`: Output directory containing the following generated files:
  - Aligned BAM files.
  - `gene_count_matrix.csv`: Gene count matrix.
  - `Covariatesstat.csv`: Metadata file containing sample information.

### 2. `pca.sh`

**Description**: This script performs Principal Component Analysis (PCA) on the gene count matrix to visualize the distribution of samples.

**Docker Input**:
- `/scratch/`: Directory containing:
  - `gene_count_matrix.csv`: Gene count matrix.
  - `Covariatesstat.csv`: Metadata file.

**Docker Output**:
- `/scratch/output/`: Output directory containing:
  - `gene_count_matrix_pca_plot.png`: PCA plot of RNA-Seq samples based on the original data.
  - `gene_count_matrix_pca_plot_filtered.png`: PCA plot based on filtered data.

### 3. `deseq2.sh`

**Description**: This script performs differential expression analysis using DESeq2 and generates results such as differential expression tables and visualizations.

**Docker Input**:
- `/scratch/`: Directory containing:
  - `gene_count_matrix.csv`: Gene count matrix.
  - `Covariatesstat.csv`: Metadata file.
  - `reference_group` and `organism`: Parameters passed via command line.

**Docker Output**:
- `/scratch/output/`: Output directory containing:
  - CSV files with differential expression results, named as `DEG_<group>_vs_<reference_group>.csv`.
  - `gene_count_matrix_filtered_count_matrix.csv`: Filtered gene count matrix.
  - `gene_count_matrix_venn_diagram.png`: Venn diagram of significant genes.

### 4. `heatmap.sh`

**Description**: This script generates heatmaps based on the filtered gene count matrix, visualizing the expression levels of significant genes across samples.

**Docker Input**:
- `/scratch/`: Directory containing:
  - `filtered_count_matrix.csv`: Filtered gene count matrix.
  - `Covariatesstat.csv`: Metadata file.

**Docker Output**:
- `/scratch/output/`: Output directory containing:
  - `gene_count_matrix_heatmap_filtered.png`: Heatmap of significant genes based on filtered data.

## Recombining Genome Files

To recombine the split genome and GTF files in the `Genome` directory, run the following command:
```bash
   ./extractAndRecombine.sh
```
This will decompress and recombine the parts into their original FASTA and GTF formats.

## Data Source

The RNA-Seq data used in this pipeline is derived from the following study:
- **Effect of deletion of CPES on gene expression during spermatogenesis in Drosophila testes**
- **SRA:** SRP387275
- **BioProject:** PRJNA860566
- **GEO:** GSE208655

### Abstract

To investigate the role of CPES in germ cell differentiation during spermatogenesis in Drosophila testis, CPES null mutants were generated and rescued with Bam-Gal4 and UAS-CPES. RNA was extracted from dissected testes and analyzed using bulk RNA-Seq on an Illumina platform. The data was processed using DESeq2 and further analyzed for gene expression changes.

## External Link

The study associated with this data is titled "Delivery of ceramide phosphoethanolamine lipids to the cleavage furrow through the endocytic pathway is essential for male meiotic cytokinesis."

## Contact

For any questions or issues with this pipeline, please open an issue on the GitHub repository or contact the repository owner.

## License

This repository is licensed under the MIT License. See the `LICENSE` file for more information.
