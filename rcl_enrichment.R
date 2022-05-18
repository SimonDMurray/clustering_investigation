# Example running
# Rscript rcl_visualisation.R --seurat=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/data/vitessce_ref.h5Seurat --rcl=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/azimuth_gex_results/rcl_mcl/rcl.res1600.txt --hallmark=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/data/hallmark.txt --output=/lustre/scratch117/cellgen/cellgeni/TIC-misc/tic-1129/test/
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library('org.Hs.eg.db')
library(SoupX)
theargs <- R.utils::commandArgs(asValues=TRUE)
input_seurat <- theargs$seurat
input_rcl <- theargs$rcl
input_hallmark <- theargs$hallmark
cutoff <- theargs$cellnumber
output_path <- theargs$output
#Reading in Seurat object and copying data matrix to count matrix
srat <- LoadH5Seurat(input_seurat)
srat@assays$RNA@counts <- srat@assays$RNA@data
#Running standard processing steps
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
srat<- FindNeighbors(srat)
#Reading in rcl table
rcl_table <- read.table(input_rcl, sep = "\t")
colnames(rcl_table) <- c("barcode", "cluster")
rcl_table <- rcl_table[rcl_table$barcode != "dummy",]
#Assigning each barcode to a rcl cluster
srat@meta.data <- cbind(srat@meta.data, rcl_table$cluster)
colnames(srat@meta.data) <- c("celltype.l1", "celltype.l2", "celltype.l3", "ori.index", "nCount_refAssay", "nFeature_refAssay", "cluster")
#Identifying clusters with more than input cell number
cluster_size <- as.data.frame(dplyr::count(srat@meta.data, cluster))
colnames(cluster_size) <- c("cluster", "cell_number")
clusters_used_df <- cluster_size[cluster_size$cell_number > int(cutoff),]
clusters_used <- as.vector(clusters_used_df$cluster)
#Subsetting seurat to cells in clusters with size > 10
barcodes_needed <- rownames(srat@meta.data[srat@meta.data$cluster %in% clusters_used,])
#Running quickmarkers
qmatrix <- srat@assays$RNA@counts
qmarkers <- quickMarkers(qmatrix, srat$cluster, N=1000, FDR=0.05)
qmarkers <- qmarkers %>% arrange(desc(tfidf))
write.table(qmarkers, paste(output_path, "qmarkers.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
#Reading in hallmark genes
gmt <- read.gmt(input_hallmark)
enrich_dict <- vector(mode="list")
enrich_dict <- vector(mode="list")
for (cluster in clusters_used) {
  cluster_markers <- qmarkers[qmarkers$cluster == cluster, ]
  cluster_genes <- as.vector(cluster_markers$gene)
  enrich_output <- enricher(cluster_genes, TERM2GENE = gmt)
  if (is.null(enrich_output)) {
    print(paste("Cluster", cluster, "does not have enough genes"))
  }
  else {
    enrich_result <- as.data.frame(enrich_output@result)
    enrich_result <- enrich_result %>% arrange(pvalue)
    enrich_dict[as.character(cluster)] <- list(enrich_result)
  }
}
#export metadata table
write.table(srat@meta.data,  paste(output_path, "metadata.txt", sep = ""), sep = "\t", quote = FALSE)
#identifying cell type and subtypes
factor_table_3_2 <- setDT(srat@meta.data)[, list(Subfactors = paste(unique(celltype.l3), collapse = " ")), by = celltype.l2]
colnames(factor_table_3_2) <- c("Type", "Subtype")
factor_table_2_1 <- setDT(srat@meta.data)[, list(Subfactors = paste(unique(celltype.l2), collapse = " ")), by = celltype.l1]
colnames(factor_table_2_1) <- c("Type", "Subtype")
#exporting subtype tables
write.table(factor_table_3_2, paste(output_path, "factor_table_3_2.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(factor_table_2_1, paste(output_path, "factor_table_2_1.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
output_r = paste(output_path, "rcl_visualisation.Rdata", sep = "")
save.image(output_r)
