library(Seurat)
library(sceasy)
library(dimRed)
theargs <- R.utils::commandArgs(asValues=TRUE)
output_path <- theargs$output
input_path <- theargs$input
do_rds <- !is.null(theargs$rds)
if (do_rds) {
  srat <- readRDS(input_path)
} else {
  adj.matrix <- Read10X(input_path)
  srat <- CreateSeuratObject(adj.matrix) 
}
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
sceasy::convertFormat(srat, from="seurat", to="anndata", outFile=paste(output_path, 'umap_seurat.h5ad', sep=""), assay = "RNA")
embed_matrix <- srat@reductions$pca@cell.embeddings
pca_df <- dimRed::as.data.frame(embed_matrix)
for(row in 1:nrow(embed_matrix)) {
  euclidean <- sqrt(sum(embed_matrix[row, ]^2))
  dispersion <- sum(embed_matrix[row, ]^2) / sum(abs(embed_matrix[row, ]))^2
  pca_df[row, "euclidean"] <- euclidean
  pca_df[row, "dispersion"] <- dispersion
}
write.table(pca_df, file=paste(output_path, "pca_with_stats.txt", sep=""), row.names=TRUE, col.names = TRUE)
