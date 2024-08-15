#!/usr/bin/env Rscript

library(DESeq2)
library(AnnotationDbi)
library(VennDiagram)
setwd("/scratch")
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript deseq2_with_venn.R <count_matrix_file> <metadata_file> <reference_group> <organism> [padj_cutoff] [log2fc_cutoff] [count_sep] [meta_sep]")
}

count_matrix_file <- args[1]
metadata_file <- args[2]
reference_group <- args[3]
organism <- args[4]
padj_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.01)
log2fc_cutoff <- ifelse(length(args) >= 6, as.numeric(args[6]), 2)
count_sep <- ifelse(length(args) >= 7, ifelse(args[7] == "tab", "\t", args[7]), ",")
meta_sep <- ifelse(length(args) >= 8, ifelse(args[8] == "tab", "\t", args[8]), ",")

count_matrix <- read.table(count_matrix_file, sep = count_sep, row.names = 1, header = TRUE)
metadata <- read.table(metadata_file, sep = meta_sep, row.names = 1, header = TRUE)

if (organism == "Homosapiens") {
  library(org.Hs.eg.db)
  ann_db <- org.Hs.eg.db
} else if (organism == "Musmusculus") {
  library(org.Mm.eg.db)
  ann_db <- org.Mm.eg.db
} else if (organism == "Drosophilamelanogaster") {
  library(org.Dm.eg.db)
  ann_db <- org.Dm.eg.db
} else {
  stop("Organismo non supportato. Usa 'Homo sapiens', 'Mus musculus' o 'Drosophila melanogaster'")
}

gene_ids <- rownames(count_matrix)
gene_symbols <- mapIds(ann_db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbols[is.na(gene_symbols)] <- "null"
annotated_gene_ids <- paste0(gene_ids, ":", gene_symbols)
rownames(count_matrix) <- annotated_gene_ids

metadata$Group <- factor(metadata$Group, levels = unique(metadata$Group))
metadata$Group <- relevel(metadata$Group, ref = reference_group)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ Group)
dds <- DESeq(dds)

results_names <- levels(metadata$Group)
contrast_groups <- setdiff(results_names, reference_group)

output_dir <- "/scratch/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Ottieni il nome base del file della matrice di conteggio
file_base_name <- tools::file_path_sans_ext(basename(count_matrix_file))

significant_genes <- list()
for (group in contrast_groups) {
  res <- results(dds, contrast = c("Group", group, reference_group))
  res <- res[order(res$padj, na.last = NA), ]
  res_filtered <- res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > log2fc_cutoff), ]
  significant_genes[[group]] <- rownames(res_filtered)
  write.csv(as.data.frame(res), file.path(output_dir, paste0(file_base_name, "_DEG_", group, "_vs_", reference_group, ".csv")))
}

filtered_count_matrix <- count_matrix[unique(unlist(significant_genes)), ]
write.table(filtered_count_matrix, file.path(output_dir, paste0(file_base_name, "_filtered_count_matrix.csv")), sep = ",", col.names = NA, quote = FALSE)

# Generazione del Venn Diagram
if (length(significant_genes) > 1) {
  venn_file <- file.path(output_dir, paste0(file_base_name, "_venn_diagram.png"))
  num_sets <- length(significant_genes)
  colors <- rainbow(num_sets)

  venn.plot <- venn.diagram(
    x = significant_genes,
    category.names = names(significant_genes),
    filename = venn_file,
    output = TRUE,
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    col = "black",
    fill = colors,
    alpha = 0.50,
    label.col = "black",
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = colors,
    cat.cex = 2.0,
    cat.fontfamily = "serif",
    margin = 0.2,
    disable.logging = TRUE
  )
}
