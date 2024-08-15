# Detect-seq: Genome-Wide Assessment of CBE Off-Targets

Detect-seq is a bioinformatics tool designed for the unbiased, genome-wide assessment of off-target effects associated with cytosine base editors (CBEs). This tool is essential for researchers in the field of gene editing, providing a comprehensive method to evaluate the specificity and safety of CRISPR-Cas9-based base editing.

## Overview

Detect-seq provides a robust pipeline for identifying and analyzing off-target mutations induced by CBEs. The workflow involves several key steps:
- Mapping sequencing reads
- Re-aligning low-quality reads
- Removing duplicates
- Detecting C-to-T conversions that signal potential off-target effects

## Input Files

- **Reference Genome**: A `.fa` or `.fasta` file containing the reference genome sequence. The script automatically detects this file in the specified genome directory.

- **FASTQ Files**: Paired-end sequencing data in `.fastq.gz` format. These files should be placed in a directory, and the path to this directory should be provided when running the script.

- **Adapters (Optional)**: Sequences for adapter trimming (default values are `AGATCGGAAGAGCACACGT` for adapter 1 and `AGATCGGAAGAGCGTCGTG` for adapter 2).

- **Threshold Value**: A numeric value used for filtering based on read quality or other criteria during processing.

## Output Files

- **BAM Files**: Aligned sequence data in BAM format, including:
  - Initial HISAT3N alignments.
  - Re-aligned BAM files after processing with BWA MEM.
  - Merged and sorted BAM files.
  - Filtered and duplicate-removed BAM files.

- **PMAT Files**: Processed mutation matrices in `.pmat` format, which contain information about detected C-to-T conversions and other mutational data.

- **Log Files**: Logs for each step of the process, which provide details on the operations performed and any potential issues encountered.

- **Filtered Output Files**: BED and WIG format files that include filtered data based on specific search strings like "CT" and "GA". These files are used for further downstream analysis.


## Usage

### Running Detect-seq

To run Detect-seq, you need to execute the following commands:

1. **Launch the Docker container with the necessary parameters**:

    ```bash
    ./run_detectseq.sh /path/to/genome_folder /path/to/fastq_folder threshold_value [adapt1] [adapt2]
    ```

    - `/path/to/genome_folder`: The directory containing the reference genome file(s) (`.fa` or `.fasta`).
    - `/path/to/fastq_folder`: The directory containing the FASTQ files.
    - `threshold_value`: The threshold value for filtering.
    - `[adapt1]` and `[adapt2]` are optional adapter sequences, with default values set to `AGATCGGAAGAGCACACGT` and `AGATCGGAAGAGCGTCGTG`, respectively.

2. **Script Workflow**:

    The script performs the following tasks:

    - Checks for existing BWA and HISAT3N indices for the reference genome.
    - Generates indices if they do not exist.
    - Processes FASTQ files by trimming adapters, aligning reads to the genome using HISAT3N, re-aligning low-quality reads with BWA MEM, and filtering the results.
    - Combines and processes the aligned data to produce BAM files and converts them to PMAT format.
    - Outputs relevant files for further analysis.

## Citation

If you use Detect-seq in your research, please cite the following paper:

Haowei Meng, Zhixin Lei, Xichen Rao, Huanan Zhao, Chengqi Yi. **Detect-seq: an unbiased method for genome-wide CBE off-targets assessment.** Nature Protocols (2023). DOI: [10.1038/s41596-023-00837-4](https://doi.org/10.1038/s41596-023-00837-4)

You should also refer to the original Detect-seq GitHub repository: [Detect-seq on GitHub](https://github.com/menghaowei/Detect-seq).

## License

This project is licensed under the Mozilla Public License Version 2.0. See the [LICENSE](https://github.com/Reproducible-Bioinformatics/CARN_Backend/blob/main/LICENSE) file for details.

## Contact

For any questions or issues, please submit an issue on the [GitHub repository](https://github.com/Reproducible-Bioinformatics/CARN_Backend).
