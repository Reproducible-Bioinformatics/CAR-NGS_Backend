# Single Cell RNA-seq Analysis Pipeline for Drosophila

This repository contains a pipeline to perform single-cell RNA sequencing (scRNA-seq) analysis for **Drosophila melanogaster** data. The raw data used in this analysis comes from the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6327103) and is organized under the `Data` directory. The pipeline is divided into multiple steps, each handled by a separate Docker container. This workflow covers the entire process from raw read alignment to downstream analysis, including differential expression and enrichment analysis.

### Project Structure

The project is organized as follows:

- **Scripts:**
  - `1_runDocker_alignIndex.sh`: Script to run the Docker container for alignment and indexing of raw sequencing data.
  - `2_runDocker_mitoRibo.sh`: Script for filtering mitochondrial and ribosomal genes.
  - `3_runDocker_clustering.sh`: Script to perform clustering analysis on the filtered data.
  - `4_runDockerFeatureselection.sh`: Script for feature selection using differential expression analysis.
  - `5_runDockerEnrichmentAnalysis.sh`: Script to run enrichment analysis on differentially expressed genes.

- **Data Directory:**
  - `Data/SRR6327103_1.fastq.gz`, `SRR6327103_2.fastq.gz`: Raw sequencing files for sample SRR6327103.
  - `Data/SRR6327104_1.fastq.gz`, `SRR6327104_2.fastq.gz`: Raw sequencing files for sample SRR6327104.
  - `combined_filtered_matrix_with_sample.mtx`, `combined_filtered_with_sample_genes.tsv`, `combined_filtered_with_sample_barcodes.tsv`: Filtered matrix, gene list, and barcodes for downstream analysis.
  - `anova_DE_results.csv`: Results from differential expression analysis.

- **Genome Directory:**
  - The genome files are divided into parts (`part_aa.gz`, `part_ab.gz`, etc.) and must be recombined using the `extractAndRecombine.sh` script. This script will generate the necessary `.fasta` and `.gtf` files for alignment.

- **Docker Information:**
  - `Docker_info/backend_alignment`: Contains the scripts for alignment and matrix combination (`index_align.sh`, `combine_matrices.R`).
  - `Docker_info/backend_DownstreamAnalysis`: Contains the R scripts for downstream analysis, including clustering (`clustering.R`), feature selection (`featureSelection.R`), and enrichment analysis (`enrichment_analysis.R`).

### Docker Containers

- **Alignment Docker** (`repbioinfo/carncellranger2`):
  This Docker container is used for aligning raw reads to the genome and indexing them for further analysis. The container runs the `index_align.sh` script and outputs the aligned BAM files.

- **Downstream Analysis Docker** (`repbioinfo/singlecelldownstream`):
  This container handles clustering, differential expression, and enrichment analysis. It runs the respective R scripts from `Docker_info/backend_DownstreamAnalysis`.

### Data Preparation

Before running the analysis, the genome files in the `Genome` directory must be extracted and recombined. Use the following command:

```bash
./extractAndRecombine.sh part_aa.gz part_ab.gz part_ac.gz part_ad.gz part_ae.gz part_af.gz
```
This will generate the `fasta` and `gtf` files necessary for alignment in the `Genome` directory. These files will be used in the alignment step, where raw sequencing reads are mapped to the Drosophila genome.

### Running the Pipeline

Each step of the pipeline is executed using a specific shell script that interacts with the appropriate Docker container. Below is a description of each script, its inputs, and outputs:

1. **Alignment and Indexing** (`1_runDocker_alignIndex.sh`):
   - **Input**: FASTQ files from the `Data` directory, genome `.fasta` and `.gtf` files from the `Genome` directory.
   - **Output**: Aligned BAM files and combined matrix of counts.
   - **Docker**: `repbioinfo/carncellranger2`.

2. **Mitochondrial and Ribosomal Gene Filtering** (`2_runDocker_mitoRibo.sh`):
   - **Input**: Aligned BAM files and matrix of counts from the previous step.
   - **Output**: Filtered matrix excluding mitochondrial and ribosomal genes.
   - **Docker**: `repbioinfo/singlecelldownstream`.

3. **Clustering** (`3_runDocker_clustering.sh`):
   - **Input**: Filtered matrix of gene expression counts.
   - **Output**: Clustering results, UMAP coordinates, and stability information.
   - **Docker**: `repbioinfo/singlecelldownstream`.

4. **Feature Selection** (`4_runDockerFeatureselection.sh`):
   - **Input**: Filtered matrix and clustering output.
   - **Output**: Differential expression results (ANOVA-like analysis).
   - **Docker**: `repbioinfo/singlecelldownstream`.

5. **Enrichment Analysis** (`5_runDockerEnrichmentAnalysis.sh`):
   - **Input**: Results from the differential expression analysis.
   - **Output**: Pathway enrichment analysis results (e.g., KEGG pathways).
   - **Docker**: `repbioinfo/carn`.

Each of these scripts is designed to be executed sequentially, with the output from one step being used as the input for the next step. The Docker containers are responsible for running the necessary bioinformatics tools and R scripts.

## Script 1: Alignment and Indexing

### Step 1: Alignment and Indexing Script

The first script handles the alignment and indexing of the FASTQ files generated from single-cell RNA-seq experiments using **Cell Ranger**. It runs within a Docker container for reproducibility and portability, ensuring all dependencies are managed.

The main steps involved in the script are:

1. **Setting up Directories**: The script defines the paths for the data (`Data`) and genome (`Genome`) directories. These directories contain the raw FASTQ files and the genome files (FASTA and GTF) required for alignment.

2. **FASTQ File Name Correction**: The script checks the names of the FASTQ files to ensure they follow the proper naming convention expected by **Cell Ranger**. If any names are inconsistent (e.g., missing lane information, sample name, or read direction), the script renames them accordingly. It logs all these changes to a log file for traceability.

3. **Genome Indexing**: The script searches for the reference genome files (FASTA and GTF) in the genome directory. If the Cell Ranger index files don't already exist, the script generates the index using the provided FASTA and GTF files. This index is stored in a directory under `Genome`.

4. **Cell Ranger Count**: Once the FASTQ files and genome index are ready, **Cell Ranger** is used to process each sample. It aligns the reads to the reference genome, generates a gene-cell count matrix, and optionally saves the BAM file if specified.

5. **Combining Matrices**: After processing each sample, the script modifies the barcodes for each cell to include the sample name and combines the count matrices from multiple samples. This step ensures that each cell is uniquely identifiable by both its barcode and sample origin.

6. **Final Output**: The combined count matrices are saved in the `Data` directory. The script also outputs modified barcodes and a log file summarizing the renaming and alignment steps.

The script operates using **Docker**, pulling a container that includes **Cell Ranger** and all required dependencies, ensuring that the environment is fully reproducible.

This script is crucial for preparing your data for downstream analysis, and the outputs include:
- Corrected FASTQ file names.
- Indexed genome reference files.
- Gene-cell count matrices for each sample.
- A combined count matrix for all samples.

### Step 2: Mitochondrial and Ribosomal Filtering Script

The second script is responsible for filtering cells based on their mitochondrial and ribosomal gene content. This is crucial for quality control in single-cell RNA-seq experiments, as cells with high mitochondrial gene expression are often stressed or dying, and those with high ribosomal gene expression might indicate contamination or technical artifacts.

#### Overview of the Script

This script runs within a **Docker** container and uses R to analyze the count matrix, identifying cells based on their mitochondrial and ribosomal gene content, then filtering out those outside specified thresholds. There are two modes of operation: handling **dense** matrices and **sparse** matrices.

The script performs the following steps:

1. **Set Up Directories**: The script sets the data directory (`Data`) containing the count matrices and sets the range of acceptable percentages for mitochondrial and ribosomal gene expression.

2. **Matrix Type Detection**: The script automatically detects whether the input matrix is **dense** (e.g., CSV) or **sparse** (Matrix Market format). For sparse matrices, it also requires the corresponding gene and barcode files.

3. **Gene Identification**: Using pattern matching, the script identifies ribosomal genes (typically prefixed with "RPS" or "RPL") and mitochondrial genes (prefixed with "MT:"). It calculates the percentage of these gene types that are expressed in each cell.

4. **Cell Filtering**: The script filters cells based on user-defined ranges for the percentage of mitochondrial and ribosomal genes expressed. Only cells within the acceptable ranges are retained in the final matrix.

5. **Plot Generation**: A scatter plot is generated to visualize the relationship between ribosomal and mitochondrial gene expression across cells. The plot helps to visualize how cells are distributed and how many fall outside the filtering thresholds.

6. **Saving Results**:
   - The filtered matrix is saved in the `Data` directory.
   - A plot of mitochondrial vs. ribosomal gene expression is saved as a PNG file.
   - The script can handle both dense and sparse matrices, saving the filtered matrix in the appropriate format (either `.mtx` for sparse or `.csv`/`.txt` for dense).

#### Inputs:
- **Dense Matrix**: A CSV file where each row represents a gene and each column represents a cell. The separator for the CSV is defined as a parameter.
- **Sparse Matrix**: A Matrix Market format file (MTX), along with the corresponding gene and barcode files.
- **Mitochondrial and Ribosomal Ranges**: The percentage ranges for filtering cells based on mitochondrial and ribosomal gene content.

#### Outputs:
- A filtered count matrix with cells that fall within the specified ranges for mitochondrial and ribosomal percentages.
- A scatter plot showing the relationship between ribosomal and mitochondrial gene expression.

This script provides a key quality control step in the single-cell RNA-seq analysis pipeline, ensuring that only high-quality cells are passed through for downstream analysis.

### Step 3: Clustering and Stability Analysis Script

The third script is used for clustering cells and assessing the stability of the clusters. It performs clustering using Seurat and calculates cluster stability via bootstrapping, ensuring that the identified clusters are consistent and reliable.

#### Overview of the Script

This script uses **Docker** to run clustering on both dense and sparse single-cell RNA-seq matrices. It also calculates cluster stability using a bootstrapping method, removing a percentage of cells in each iteration and recalculating the clusters. Stability scores are derived by comparing the clusters before and after bootstrapping using the **Jaccard Index**.

Key steps of the script include:

1. **Set Up Directories**: The script points to the `Data` directory, where the count matrices are stored, and defines the clustering parameters such as bootstrap percentage, stability threshold, number of permutations, and Seurat clustering resolution.

2. **Matrix Type Detection**: The script detects whether the input matrix is **dense** (CSV format) or **sparse** (Matrix Market format). For sparse matrices, it also requires the corresponding gene and barcode files.

3. **Custom Clustering Function**: The script applies a custom clustering function using **Seurat**, a widely used R package for single-cell analysis. This function creates a Seurat object, normalizes the data, identifies variable features, scales the data, runs PCA, finds neighbors, and performs clustering. It also performs dimensionality reduction using **UMAP** to visualize the clusters.

4. **Stability Analysis**: The script performs stability analysis using bootstrapping. In each iteration, a percentage of cells are removed from the dataset, and clustering is recalculated. The **Jaccard Index** is used to compare the clusters before and after bootstrapping. A stability score is assigned to each cell based on how often it remains in the same cluster across iterations.

5. **Plot Generation**: Two plots are generated:
   - A **UMAP plot** colored by cluster.
   - A **UMAP plot** colored by stability scores, with a fixed range from 0 to 1.

6. **Saving Results**:
   - The clustering results, including the cluster assignments and UMAP coordinates, are saved as a CSV file.
   - Stability scores are included in the output, providing a measure of the consistency of each cell's cluster.
   - UMAP plots are saved as PNG files to visually assess cluster assignments and stability.

#### Inputs:
- **Dense Matrix**: A CSV file where rows represent genes and columns represent cells. The separator for the CSV is defined as a parameter.
- **Sparse Matrix**: A Matrix Market format file (MTX) with the corresponding gene and barcode files.
- **Bootstrap Percentage**: The percentage of cells to remove in each bootstrap iteration (e.g., 10%).
- **Stability Threshold**: The minimum Jaccard Index value for a cluster to be considered stable (e.g., 0.8).
- **Permutations**: The number of bootstrap iterations to perform (e.g., 10).
- **Seurat Resolution**: The resolution parameter for Seurat clustering, which controls the granularity of the clusters (e.g., 0.8).

#### Outputs:
- A CSV file containing cluster assignments, UMAP coordinates, and stability scores for each cell.
- Two UMAP plots:
  - One colored by cluster.
  - One colored by stability scores, allowing for a quick assessment of the stability of the clusters.
- The filtered matrix is saved in the `Data` directory.

This script plays a critical role in identifying stable, reliable clusters from single-cell RNA-seq data and visualizing the results for further interpretation.
### Step 4: Feature Selection Script

The fourth script focuses on **feature selection**, specifically identifying differentially expressed (DE) genes using various statistical methods like **ANOVA**, **MAST**, and **edgeR**. This script uses clustering information from the previous step and filters genes based on thresholds such as log2 fold change (log2FC) and p-value, which are commonly used to identify genes with significant differences between clusters.

#### Overview of the Script

This script identifies DE genes across cell clusters. It can handle both **dense** (CSV format) and **sparse** (Matrix Market format) input matrices. It also generates heatmaps and volcano plots to visualize the DE genes.

Key steps of the script include:

1. **Set Up Directories**: The script points to the `Data` directory where the count matrices and clustering files are stored. It defines several parameters for feature selection, including log2FC and p-value thresholds.

2. **Matrix Type Detection**: The script detects whether the input matrix is **dense** (CSV format) or **sparse** (Matrix Market format). For sparse matrices, it requires the corresponding gene and barcode files.

3. **ANOVA and MAST**: The script performs **ANOVA** using **MAST**, which is designed for single-cell RNA-seq data, and **edgeR**, a popular package for RNA-seq differential expression analysis. The results include log2 fold change, p-values, and DE gene identification.

4. **Volcano Plot**: A volcano plot is generated to visualize the DE genes, plotting log2FC on the x-axis and -log10(p-value) on the y-axis. This plot highlights genes that are significantly upregulated or downregulated in each cluster comparison.

5. **Heatmap Generation**: The script generates heatmaps for DE genes. This helps in visualizing how these genes are expressed across the different cells and clusters.

6. **Pairwise Cluster Comparisons**: The script performs pairwise comparisons between clusters, identifying genes that are differentially expressed between each pair of clusters.

7. **Stability Filtering**: The script can filter cells based on the stability scores calculated in the previous step. Only stable cells (those with a stability score greater than a defined threshold) are included in the feature selection analysis.

#### Inputs:
- **Dense Matrix**: A CSV file where rows represent genes and columns represent cells. The separator for the CSV is defined as a parameter.
- **Sparse Matrix**: A Matrix Market format file (MTX) with the corresponding gene and barcode files.
- **Clustering File**: A CSV file containing the clustering results (with columns for `Cell`, `Cluster`, and `Stability`).
- **Log2FC Threshold**: The log2 fold change threshold for identifying DE genes (e.g., 1).
- **P-value Threshold**: The p-value threshold for identifying DE genes (e.g., 0.05).
- **Stability Threshold**: The stability threshold for filtering cells based on their stability score (e.g., 0.8).

#### Outputs:
- A CSV file containing DE gene results for each comparison.
- **Volcano plots** for each comparison, showing upregulated and downregulated genes.
- **Heatmaps** for visualizing DE gene expression across clusters.
- The script can be set to perform comparisons for all cells, as well as comparisons that include only **stable cells** (those passing the stability threshold).

This script is critical for identifying key genes that differentiate cell populations in single-cell RNA-seq data. The visualizations (volcano plots and heatmaps) provide an intuitive way to assess which genes are driving the differences between cell clusters.

### Step 5: Enrichment Analysis Script

The final step in the pipeline involves running **enrichment analysis** on the differentially expressed (DE) genes obtained from previous steps. This script identifies biological pathways or processes in which the DE genes are involved using sources like **KEGG**, **GO**, or **Reactome**.

#### Overview of the Script

The enrichment analysis script processes the results of differential expression and performs pathway enrichment analysis. This analysis identifies overrepresented biological processes or pathways in the DE genes using external databases.

Key steps of the script include:

1. **Data Preparation**:
    - The script begins by setting up the working environment and specifying the DE gene file (`anova_DE_results.csv`), species (`dmelanogaster` for *Drosophila melanogaster*), and the source of enrichment (e.g., **KEGG** or **GO**).
    - It ensures that the DE gene file contains the necessary columns (`logFC`, `PValue`, and `FDR`) and prepares gene names for analysis, converting them to gene symbols if needed.

2. **Enrichment Analysis**:
    - The script runs an enrichment analysis using the **g:Profiler2** package. It sends the list of DE genes to external databases like **KEGG** or **GO** to find overrepresented pathways or biological processes.
    - A **bar plot** is generated to visualize the most enriched terms based on fold enrichment and the significance of the terms (-log10(p-value)).

3. **Dynamic Plot Generation**:
    - The plot visualizes up to 20 enriched terms and presents the **fold enrichment** on the x-axis and the **-log10(p-value)** on the y-axis, providing an intuitive way to assess the significance of each pathway.
    - The plot is saved as a **PDF** file for easy sharing and publication.

4. **Input Parameters**:
    - **DE File**: This is a CSV file containing the results of differential expression analysis. It includes gene names, log fold change, p-values, and FDR values.
    - **Species**: The species is defined as `dmelanogaster` for *Drosophila melanogaster*. The script can be adapted to other species by modifying the species argument.
    - **Source**: The source of enrichment data can be **KEGG**, **GO** (Biological Process, Molecular Function, Cellular Component), **Reactome**, or **WikiPathways**.
    - **Separator**: The separator in the DE gene file, such as a comma (`,`), tab (`\t`), etc.
    - **Max Terms**: The maximum number of enriched terms to display in the output plot (e.g., 20 terms).

5. **Outputs**:
    - **Enrichment Results**: A PDF file containing a bar plot of the enriched terms.
    - **Error Handling**: The script provides error messages if no significant genes are found or if the enrichment analysis does not return any results.

This script completes the single-cell RNA-seq analysis by contextualizing DE genes in the larger framework of biological pathways and processes. The visualizations provide insights into the key pathways active in different cell clusters.
