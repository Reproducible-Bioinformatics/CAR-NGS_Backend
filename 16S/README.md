# 16S rRNA Bioinformatics Analysis Pipeline

This repository contains a pipeline for the analysis of 16S rRNA gene sequencing data, typically used for microbial community profiling. The pipeline is designed to process paired-end FASTQ files, perform quality control, denoise the sequences using DADA2, assign taxonomy using the Silva reference database, and generate various outputs such as feature tables, taxonomy summaries, and visualizations.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Overview

The 16S rRNA Bioinformatics Analysis Pipeline processes raw sequencing data from 16S rRNA gene amplicons. The pipeline includes the following key steps:

1. **Quality Control**: Checks the quality of the sequencing reads and removes low-quality bases.
2. **Denoising**: Uses DADA2 to correct errors in sequencing data and infer true biological sequences.
3. **Taxonomic Classification**: Assigns taxonomy to the inferred sequences using a pre-trained Silva 138 classifier.
4. **Visualization**: Generates visualizations of the taxonomic composition of the samples.

This pipeline is designed to be run within a Docker container to ensure reproducibility and ease of use.

## Requirements

- Docker installed on your system
- Paired-end FASTQ files from 16S rRNA sequencing
- Pre-trained Silva 138 classifier in QIIME2 format (`.qza`)

## Installation

1. Clone this repository to your local machine:
    ```bash
    git clone https://github.com/yReproducible-Bioinformatics/CAR-NGS_Backend
    cd 16S-pipeline
    ```

2. Ensure that you have Docker installed and running on your system.

## Usage

To run the pipeline, place your paired-end FASTQ files in a directory and specify this directory as `LOCAL_FASTQ_DIR` when running the Docker container.

### Running the Pipeline

```bash
docker run --rm -it \
-v ${LOCAL_FASTQ_DIR}:/scratch \
--user $(id -u):$(id -g) \
repbioinfo/qiime2023 /home/qiime_full.sh
```
### Example

```bash
docker run --rm -it \
-v /path/to/your/fastq:/scratch \
--user $(id -u):$(id -g) \
repbioinfo/qiime2023 /home/qiime_full.sh

```
This command will process all paired-end FASTQ files in the specified directory and output the results to the /scratch/aligned_results directory within the container.
## Pipeline Steps

1. **Data Import**: The pipeline imports paired-end sequences into a QIIME2 artifact format (`.qza`).
2. **Demultiplexing Summary**: A summary of the sequencing data is generated to assess quality.
3. **Denoising with DADA2**: The data is denoised to remove errors and infer ASVs (Amplicon Sequence Variants).
4. **Taxonomic Classification**: Taxonomy is assigned to the ASVs using the Silva 138 reference database.
5. **Visualization**: The taxonomic composition is visualized using bar plots.

## Output Files

The pipeline generates the following key output files:

- **paired-end-demux.qza**: QIIME2 artifact containing the demultiplexed sequences.
- **paired-end-demux.qzv**: Visualization of the demultiplexed sequences.
- **table.qza**: Feature table artifact containing the abundance of ASVs.
- **rep-seqs.qza**: Representative sequences of ASVs.
- **taxonomy.qza**: Taxonomic classifications of the representative sequences.
- **taxa-bar-plots.qzv**: Visualization of taxonomic composition as bar plots.
- **HTML Reports**: Exported visualizations and reports in HTML format for easy sharing.

## Troubleshooting

- Ensure that your FASTQ files are named with the `_R1` and `_R2` suffixes for paired-end reads.
- Verify that the path to the Silva classifier is correct and accessible within the Docker container.
- If you encounter permission issues, make sure to use the `--user $(id -u):$(id -g)` flag when running Docker.

## Contributing

Contributions are welcome! Please fork this repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

