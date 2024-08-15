#!/usr/bin/env Rscript

# Carica le librerie necessarie
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(AnnotationDbi)
library(VennDiagram)

# Leggi gli argomenti dalla linea di comando
args <- commandArgs(trailingOnly = TRUE)

# Controlla se sono stati passati i file necessari
if (length(args) < 4) {
  stop("Usage: Rscript script.R <count_matrix_file> <metadata_file> <reference_group> <organism> [count_sep] [meta_sep]")
}

# Parametri di input
count_matrix_file <- file.path("/scratch/", args[1])
metadata_file <- file.path("/scratch/", args[2])
reference_group <- args[3]
organism <- args[4]
padj_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.01)
log2fc_cutoff <- ifelse(length(args) >= 6, as.numeric(args[6]), 2)
count_sep <- ifelse(length(args) >= 7, ifelse(args[7] == "tab", "\t", args[7]), ",")
meta_sep <- ifelse(length(args) >= 8, ifelse(args[8] == "tab", "\t", args[8]), ",")

# Carica la matrice di conteggio
count_matrix <- read.table(count_matrix_file, sep = count_sep, row.names = 1, header = TRUE)

# Carica il file di metadata
metadata <- read.table(metadata_file, sep = meta_sep, row.names = 1, header = TRUE)

# Verifica che i campioni nel metadata corrispondano ai campioni nella matrice di conteggio
if (!all(rownames(metadata) %in% colnames(count_matrix))) {
  stop("I campioni nel file di metadata non corrispondono ai campioni nella matrice di conteggio")
}

# Scegli la libreria di annotazione in base all'organismo
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

# Annotazione dei geni
gene_ids <- rownames(count_matrix)
gene_symbols <- mapIds(ann_db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Sostituisci i NA con "null"
gene_symbols[is.na(gene_symbols)] <- "null"

# Crea una nuova matrice con GENEID:geneNAME
annotated_gene_ids <- paste0(gene_ids, ":", gene_symbols)
rownames(count_matrix) <- annotated_gene_ids

# Assicurati che il gruppo di riferimento esista nel metadata
if (!reference_group %in% metadata$Group) {
  stop(paste("Il gruppo di riferimento", reference_group, "non esiste nel metadata"))
}

# Esegui la PCA
pca_result <- prcomp(t(log2(count_matrix + 1)))

# Prepara i dati per il plot
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)
metadata$Sample <- rownames(metadata)
pca_data <- merge(pca_data, metadata, by = "Sample")

# Crea la cartella di output se non esiste
output_dir <- "/scratch/output/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Crea il plot PCA
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of RNA-Seq Data",
       x = paste0("PC1: ", round(100 * summary(pca_result)$importance[2, 1], 2), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca_result)$importance[2, 2], 2), "% variance")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# Salva il plot PCA nella cartella di output
ggsave(filename = file.path(output_dir, "pca_plot.png"), plot = p, width = 8, height = 6)

# Analisi dell'espressione differenziale con DESeq2
# Assicurati che l'ordine dei campioni nel metadata corrisponda a quello della matrice di conteggio
metadata <- metadata[match(colnames(count_matrix), rownames(metadata)), ]

# Converti la colonna Group in un fattore e imposta il livello di riferimento
metadata$Group <- factor(metadata$Group, levels = unique(metadata$Group))
metadata$Group <- relevel(metadata$Group, ref = reference_group)

# Crea un oggetto DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ Group)

# Fai partire l'analisi DESeq2
dds <- DESeq(dds)

# Ottieni i confronti
results_names <- levels(metadata$Group)
contrast_groups <- setdiff(results_names, reference_group)
significant_genes <- list()

# Esegui confronti con il gruppo di riferimento e salva i risultati
for (group in contrast_groups) {
  res <- results(dds, contrast = c("Group", group, reference_group))
  res <- res[order(res$padj, na.last = NA), ]
  res_filtered <- res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > log2fc_cutoff), ]

  # Aggiungi i geni significativi alla lista
  significant_genes[[group]] <- rownames(res_filtered)

  print(res)
  write.csv(as.data.frame(res), file.path(output_dir, paste0("DEG_", group, "_vs_", reference_group, ".csv")))
}

# Filtra la matrice di conteggio con i geni significativi
filtered_count_matrix <- count_matrix[unique(unlist(significant_genes)), ]
write.table(filtered_count_matrix, file.path(output_dir, "filtered_count_matrix.csv"), sep = ",", col.names = NA, quote = FALSE)

# Esegui PCA sui dati filtrati
pca_result_filtered <- prcomp(t(log2(filtered_count_matrix + 1)))
pca_data_filtered <- as.data.frame(pca_result_filtered$x)
pca_data_filtered$Sample <- rownames(pca_data_filtered)
pca_data_filtered <- merge(pca_data_filtered, metadata, by = "Sample")

# Crea il plot PCA sui dati filtrati
p_filtered <- ggplot(pca_data_filtered, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of RNA-Seq Data (Filtered by DEG)",
       x = paste0("PC1: ", round(100 * summary(pca_result_filtered)$importance[2, 1], 2), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca_result_filtered)$importance[2, 2], 2), "% variance")) +
  theme_minimal() +
  theme(legend.title = element_blank())
ggsave(filename = file.path(output_dir, "pca_plot_filtered.png"), plot = p_filtered, width = 8, height = 6)
if(length(which(rowSums(filtered_count_matrix)==0))!=0){
filtered_count_matrix=filtered_count_matrix[-which(rowSums(filtered_count_matrix)==0),]
}
# Genera e salva la heatmap con i nomi dei geni
heatmap_file <- file.path(output_dir, "heatmap_filtered.png")
pheatmap(
  mat = filtered_count_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = metadata[, "Group", drop = FALSE],
  main = "Heatmap of Significant Genes",
  filename = heatmap_file
)

# Genera e salva i Venn diagram
if (length(significant_genes) > 1) {
  venn_file <- file.path(output_dir, "venn_diagram.png")

  # Imposta un numero di colori pari al numero di set nel diagramma Venn
  num_sets <- length(significant_genes)
  #colors <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3")[1:num_sets]
colors=rainbow(num_sets)
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
    disable.logging=TRUE,
  )
}
