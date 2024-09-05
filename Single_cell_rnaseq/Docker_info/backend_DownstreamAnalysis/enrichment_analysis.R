# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gprofiler2)
library(reshape2)
library(biomaRt)
setwd("/scratch")
# Function to dynamically load the correct org.db package based on species
load_org_db <- function(species) {
  species_db <- switch(species,
                       "hsapiens" = "org.Hs.eg.db",
                       "mmusculus" = "org.Mm.eg.db",
                       "dmelanogaster" = "org.Dm.eg.db",
                       "rnorvegicus" = "org.Rn.eg.db",
                       "scerevisiae" = "org.Sc.sgd.db",
                       stop("Species not supported. Supported species: hsapiens, mmusculus, dmelanogaster, rnorvegicus, scerevisiae"))

  if (!requireNamespace(species_db, quietly = TRUE)) {
    BiocManager::install(species_db)
  }
  library(species_db, character.only = TRUE)
}

# Function to create a unique output directory
create_output_dir <- function(base_name) {
  dir_name <- paste0(base_name, "_enrichment_analysis")
  i <- 1
  final_dir <- dir_name

  while (dir.exists(final_dir)) {
    final_dir <- paste0(dir_name, "_", i)
    i <- i + 1
  }

  dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
  return(final_dir)
}

# Function to load and check data
load_and_check_data <- function(file_path, separator) {
  if (separator == "tab") {
    separator <- "\t"
  }

  data <- read.table(file_path, header = TRUE, sep = separator)

  # Check if the file has the required columns
  required_columns <- c("logFC", "PValue", "FDR")
  if (!all(required_columns %in% colnames(data))) {
    stop("The DE file must contain the columns: logFC, PValue, FDR.")
  }

  # Preliminary check on gene names
  gene_column <- data[, 1]
  prefixes <- substr(gene_column, 1, 2)

  # If gene names contain ":", extract the part after ":"
  if (all(grepl(":", gene_column))) {
    print("Annotation structure: ID:geneName")
    print("Performing change in data structure to geneName only")
    gene_column <- sapply(gene_column, FUN = function(x) strsplit(x, ":")[[1]][2])  # Extract geneName
    data[, 1] <- gene_column
  } else if (length(unique(prefixes)) <= 10) {
    print("Probably you have gene IDs. Trying to convert them to gene names.")
    mart <- useEnsembl(biomart = "ensembl", dataset = paste0(species, "_gene_ensembl"))
    conversion <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"), values = gene_column, mart = mart)
    if (nrow(conversion) > 0) {
      data <- merge(data, conversion, by.x = colnames(data)[1], by.y = "ensembl_gene_id", all.x = TRUE)
      data[, 1] <- data$external_gene_name  # Update with external_gene_name
      data <- data[, !(colnames(data) %in% "external_gene_name")]
    } else {
      print("No matching gene names found for IDs.")
    }
  } else {
    print("Assuming gene names are already correct.")
  }

  return(data)
}

# Function to create enrichment plot with dynamic title and legend, and exclude NA values
create_enrichment_plot <- function(gp_df, source, output_dir, max_terms) {
  gp_df <- gp_df[1:max_terms, c("term_name", "intersection_size", "p_value")]
  colnames(gp_df) <- c("biological_process", "fold_enrichment", "log_p_value")
  gp_df$log_p_value <- -log10(gp_df$log_p_value)  # Convert p_value to log10(p_value)

  # Remove rows with NA values
  gp_df <- na.omit(gp_df)

  enrichment_title <- switch(source,
                             "GO:BP" = "GO Biological Process Enrichment",
                             "GO:MF" = "GO Molecular Function Enrichment",
                             "GO:CC" = "GO Cellular Component Enrichment",
                             "KEGG" = "KEGG Pathway Enrichment",
                             "REAC" = "Reactome Pathway Enrichment",
                             "WP" = "WikiPathways Enrichment",
                             paste(source, "Enrichment", sep = " "))

  enrichment_plot <- ggplot(gp_df, aes(x = reorder(biological_process, fold_enrichment), y = fold_enrichment, fill = log_p_value)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "white", high = "red") +  # Set gradient from white to red
    labs(title = enrichment_title, x = "Biological Process", y = "Fold Enrichment", fill = "-log10(p-value)") +
    theme_minimal() +
    coord_flip()

  # Save the enrichment plot as a PDF
  output_file <- file.path(output_dir, "enrichment_plot.pdf")
  ggsave(filename = output_file, plot = enrichment_plot, width = 10, height = 8)
}

# Function to run enrichment analysis
run_enrichment_analysis <- function(genes, species, source, output_dir, max_terms) {
  tryCatch({
    gp_result <- gost(query = genes, organism = species, sources = source)

    # Check if there are results
    if (is.null(gp_result$result) || nrow(gp_result$result) == 0) {
      stop("No enrichment results found.")
    } else {
      create_enrichment_plot(gp_result$result, source, output_dir, max_terms)
    }
  }, error = function(e) {
    cat("Error during enrichment analysis:", e$message, "\n")
  })
}

# Main function to run the entire enrichment analysis workflow
perform_enrichment <- function(file_path, species, source, separator, max_terms) {
  # Create a unique output directory
  base_name <- tools::file_path_sans_ext(basename(file_path))
  output_dir <- create_output_dir(base_name)

  # Load and validate the data
  data <- load_and_check_data(file_path, separator)

  # Extract significant genes (e.g., FDR < 0.05)
  significant_genes <- data[data$FDR < 0.05, 1]

  # Check if there are significant genes
  if (length(significant_genes) == 0) {
    stop("No significant genes found. Check your FDR threshold.")
  }

  # Run the enrichment analysis
  run_enrichment_analysis(significant_genes, species, source, output_dir, max_terms)

  cat("Analysis completed. Results saved in:", output_dir, "\n")
}

# Call the function from command line
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]  # Path to the DE file
species <- args[2]    # Species, e.g., "hsapiens", "mmusculus"
source <- args[3]     # Source, e.g., "GO:BP", "KEGG"
separator <- args[4]  # Separator, e.g., "tab", ","
max_terms <- as.numeric(args[5])  # Maximum number of terms to show in the plot
file_path=paste("/scratch/",file_path,sep="")
perform_enrichment(file_path, species, source, separator, max_terms)
