# ATAC-seq Analysis Pipeline

This repository provides a Dockerized pipeline for performing ATAC-seq (Assay for Transposase-Accessible Chromatin with high-throughput sequencing) analysis. The pipeline automates quality control, alignment, peak calling, and conversion of results into BigWig format.

## Prerequisites

- **Docker**: Ensure Docker is installed on your system. You can download it from [Docker's official website](https://www.docker.com/products/docker-desktop).
- **FASTQ files**: Paired-end or single-end FASTQ files for ATAC-seq analysis.
- **Reference genome**: A FASTA file of the reference genome to be used for alignment.

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/atacseq-pipeline.git
cd atacseq-pipeline
```
### 2. Prepare Your Data

- Place your **FASTQ files** in a directory, for example, `/path/to/fastq_files`. These files should be named appropriately, such as `sample_R1.fastq.gz` and `sample_R2.fastq.gz` for paired-end data, or `sample.fastq.gz` for single-end data.
- Place your **reference genome FASTA file** in another directory, for example, `/path/to/fasta_files`. This directory should contain the `.fa` or `.fasta` file and any associated index files if available.

### 3. Run the Analysis

To execute the ATAC-seq analysis pipeline, use the following command:

```bash
./atacseq-pipeline.sh /path/to/fastq_files /path/to/fasta_files [THREADS]
```
- `/path/to/fastq_files`: Path to the directory containing your FASTQ files. These files should be named appropriately, such as `sample_R1.fastq.gz` and `sample_R2.fastq.gz` for paired-end data, or `sample.fastq.gz` for single-end data.
- `/path/to/fasta_files`: Path to the directory containing the FASTA file of the reference genome. This directory should contain the `.fa` or `.fasta` file and any associated index files if available.
- `[THREADS]`: (Optional) Number of threads to use for the analysis. Defaults to `8` if not provided.

### Example

```bash
./atacseq-pipeline.sh /data/fastq /data/genomes 8
```
This command will mount the specified directories into the Docker container, execute the ATAC-seq analysis, and save the results in /path/to/fastq_files/results.

## Output

The results of the analysis will be stored in the `results` subdirectory within your FASTQ files directory. The output includes:

- **Quality Control Reports**: FastQC reports for each FASTQ file, providing insights into the quality of your sequencing data.
- **BAM Files**: Aligned reads in BAM format. There will be two versions:
  - The original sorted BAM file, which includes all aligned reads.
  - A filtered BAM file (`sorted.noorg.bam`), where organellar reads (such as mitochondrial and chloroplast DNA) have been removed.
- **Index Files**: Index files (`.bai`) corresponding to each BAM file for efficient data access.
- **Peak Calling Results**: Output from MACS2, including:
  - Peak files in narrowPeak or broadPeak format.
  - BedGraph files for visualizing signal intensity across the genome.
- **BigWig Files**: Processed BedGraph files converted into BigWig format for efficient visualization in genome browsers.

## Troubleshooting

- **Docker Permission Issues**: If you encounter permission errors while running Docker, try using `sudo` before your commands.
- **Input File Issues**: Ensure that your FASTQ and FASTA files are correctly formatted and placed in the specified directories. Misnamed files or incorrect paths can lead to errors.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Acknowledgements

This pipeline leverages several bioinformatics tools, including:

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for sequence alignment.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control.
- [MACS2](https://github.com/macs3-project/MACS) for peak calling.

## Contributing

Contributions are welcome! If you'd like to contribute, please fork the repository, create a new branch, and submit a Pull Request. Issues can be reported via the GitHub Issues page.

## Contact

For questions, feedback, or issues, please open an issue on this repository or contact the maintainer at [your.email@example.com](mailto:your.email@example.com).

## Cite
GSM6929262: Dmel_78-Pupa_bioRep-1; Drosophila melanogaster; ATAC-seq (SRR23061259)
