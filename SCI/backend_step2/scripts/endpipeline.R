#!/usr/bin/env Rscript

#demultiplexing("docker", data.folder="/20tb/ratto/sciformatrix006/", threads = 1)

# This is the end of the pipeline.
# The data is now ready to load into R as a Monocle CellDataSet.
#
# Here is an R function to do that.
# Now would also be a good time to check that you are using the latest version of Monocle.
# As of December 2017, that is Monocle version 2.6.1.
#
# Also make sure to get irlba version 2.3.2 from Github (version 2.3.1 is buggy).
# devtools::install_github("bwlewis/irlba", ref = "9dab7ed2152c42e5c99a3b01e1ac33ba47e3a909")
#
#install.packages("pkgload")
library(pkgload)
#devtools::install_github("bwlewis/irlba")
#install.packages("Matrix")
library(Matrix)
#install.packages("Biobase")
library(Biobase)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("monocle")
library(monocle)

load.cds = function(mat.path, gene.annotation.path, cell.annotation.path) {
  df = read.table(
    mat.path,
    col.names = c("gene.idx", "cell.idx", "count"),
    colClasses = c("integer", "integer", "integer"))
  
  #df$gene.idx[df$gene.idx == 0] <- 1
  
  gene.annotations = read.table(
    gene.annotation.path,
    col.names = c("id", "gene_short_name"),
    colClasses = c("character", "character"))
  
  cell.annotations = read.table(
    cell.annotation.path,
    col.names = c("cell", "sample"),
    colClasses = c("character", "factor"))
  
  rownames(gene.annotations) = gene.annotations$id
  rownames(cell.annotations) = cell.annotations$cell
  
  # add a dummy cell to ensure that all genes are included in the matrix
  # even if a gene isn't expressed in any cell
  df = rbind(df, data.frame(
    gene.idx = c(1, nrow(gene.annotations)),
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    count = c(1, 1)))
  
  mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
  mat = mat[, 1:(ncol(mat)-1)]
  
  rownames(mat) = gene.annotations$id
  colnames(mat) = cell.annotations$cell
  
  pd = new("AnnotatedDataFrame", data = cell.annotations)
  fd = new("AnnotatedDataFrame", data = gene.annotations)
  
  cds = newCellDataSet(mat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
  pData(cds)$n.umi = apply(exprs(cds), 2, sum)
  
  return(cds)
}


#folder = "/20tb/ratto/Bertero/sci/sciformatrix004/final-output/"
folder = "/data/scratch/final-output/"
mat = load.cds(mat.path = paste(folder, "UMI.count.matrix", sep = ""),
               gene.annotation.path = paste(folder, "gene.annotations", sep = ""),
               cell.annotation.path = paste(folder, "cell.annotations", sep = ""))
df = as.data.frame(as.matrix(as.matrix(mat)))
df_filtered <- df[rowSums(df != 0) > 0, ]
write.table(df, paste(folder, "exp_mat.csv", sep = ""), sep = ",", row.names = TRUE, col.names = TRUE)
saveRDS(df, file = paste(folder, "exp_mat.rds", sep = ""))
write.table(df_filtered, paste(folder, "exp_mat_no0.csv", sep = ""), sep = ",", row.names = TRUE, col.names = TRUE)
saveRDS(df_filtered, file = paste(folder, "exp_mat_no0.rds", sep = ""))

# df2 = read.csv(paste(folder, "exp_mat.csv", sep = ""))
# 
# mat.path = paste(folder, "UMI.count.matrix", sep = "")
# gene.annotation.path = paste(folder, "gene.annotations", sep = "")
# cell.annotation.path = paste(folder, "cell.annotations", sep = "")

# silencing_matrix <- read.csv("/20tb/ratto/Bertero/sci/new_8_sc_5_UMI/silencing_matrix.csv", header = T, row.names = 1)
# genes = as.data.frame(row.names(silencing_matrix))
# sum(duplicated(genes$`row.names(silencing_matrix)`))
# new_row_names <- gsub("\\..*", "", rownames(silencing_matrix))
# patterns = new_row_names[duplicated(new_row_names)]
# df_filtered <- silencing_matrix[grepl(paste(patterns, collapse = "|"), rownames(silencing_matrix)), ]


