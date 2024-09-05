# **C**ross-platform **A**ccessible **R**eproducible-NGS data analysis (**CAR-NGS**)

<img src="./Logo.png" alt="Logo" width="200"/>

### "CAR-NGS: Buckle Up, We're Taking Your DNA for a Ride!"

This repository contains a collection of automated data analysis pipelines for various genomic techniques. Each pipeline is designed to simplify the process of analyzing high-throughput sequencing data using Docker, ensuring reproducibility and ease of use.

## Repository Structure

The repository is organized into directories, each corresponding to a specific genomic technique. Below is a brief description of each pipeline:

### 1. `16S`

This pipeline is designed for the analysis of 16S rRNA gene sequencing data, typically used for microbial community profiling. The pipeline includes steps for quality control, taxonomic classification, and diversity analysis.

### 2. `Atac-seq`

The ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) pipeline facilitates the identification of open chromatin regions across the genome. The pipeline includes steps for quality control, alignment, peak calling, and visualization.

### 3. `Bulk_rnaseq`

This pipeline processes bulk RNA-Seq data to analyze gene expression levels in a population of cells. The pipeline includes steps for read alignment, quantification, differential expression analysis, and data visualization.

### 4. `Detect-seq`

Detect-seq is a pipeline for detecting sequence variants from high-throughput sequencing data. The pipeline includes steps for read alignment, variant calling, and annotation.

### 5. `Gro-seq`

The GRO-seq (Global Run-On Sequencing) pipeline analyzes nascent RNA transcripts to provide insights into transcriptional activity. The pipeline includes steps for alignment, transcript assembly, and expression quantification.

### 6. `HTGTS`

High-Throughput Genome-Wide Translocation Sequencing (HTGTS) is a technique used to map chromosomal translocations. This pipeline includes steps for alignment, translocation detection, and visualization.

### 7. `SCI`

The SCI (Single Cell Indexing) pipeline processes single-cell sequencing data to analyze gene expression at the individual cell level. The pipeline includes steps for cell demultiplexing, alignment, gene quantification, and clustering analysis.

### 8. `Single_cell_rnaseq`

This pipeline is tailored for single-cell RNA-Seq data, enabling the analysis of gene expression at the single-cell level. The pipeline includes steps for quality control, normalization, clustering, and differential expression analysis.

### 9. `SpatialTranscriptomics`

The Spatial Transcriptomics pipeline analyzes data from spatially resolved transcriptomics experiments, which retain the spatial context of gene expression. The pipeline includes steps for data alignment, spatial mapping, and visualization.

### 10. `SRA_toolkit`

The SRA Toolkit pipeline provides tools for accessing and processing sequencing data stored in the Sequence Read Archive (SRA). This includes tools for downloading, converting, and processing SRA files.

### 11. `TCR_SingleCell_RNAseq`

This pipeline is focused on T-cell receptor (TCR) sequencing data from single-cell RNA-Seq experiments, which helps in understanding TCR diversity and clonal expansion. The pipeline includes steps for TCR sequence identification, annotation, and analysis.

### 12. `WholeGenomeSequencing`

The Whole Genome Sequencing (WGS) pipeline processes data from WGS experiments, enabling comprehensive analysis of genomic variants. The pipeline includes steps for alignment, variant calling, annotation, and visualization.

## How to Use These Pipelines

Each directory contains a dedicated pipeline for a specific genomic analysis technique. To run a pipeline:

1. **Navigate to the appropriate directory**:
```bash
  cd Bulk_rnaseq
```

2. **Read the specific `README.md`**: Each directory contains its own `README.md` file with detailed instructions on how to run the pipeline, including dependencies, input data formats, and output file descriptions.

3. **Prepare your environment**: Ensure Docker is installed on your system. The pipelines are designed to run inside Docker containers, ensuring consistency and reproducibility across different computing environments.

4. **Execute the pipeline**: Follow the instructions in the respective directory to run the analysis pipeline using Docker. Most pipelines can be started with a simple command, such as:
```bash
./run_pipeline.sh
```

## Prerequisites

- **Docker**: Ensure that Docker is installed and running on your system. Docker containers are used to encapsulate the environment needed for each pipeline, making the setup process easier and ensuring reproducibility.
- **Sequencing Data**: Input data formats vary depending on the pipeline, but generally include FASTQ files, BAM files, or VCF files.

## Contributing

Contributions to this repository are welcome. If you have improvements, bug fixes, or new pipelines to add, please submit a pull request. Ensure that all pipelines are well-documented and tested.

## License

This repository is licensed under the MIT License. See the `LICENSE` file for more information.

## Contact

For any questions or issues related to these pipelines, please open an issue on the GitHub repository or contact the repository owner.
