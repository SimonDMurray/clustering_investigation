library(Seurat)
library(Matrix)
theargs <- R.utils::commandArgs(asValues=TRUE)
do_neighbours <- !is.null(theargs$neighbours)
do_rds <- !is.null(theargs$rds)
output_path <- theargs$output
input_path <- theargs$input
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
if (do_neighbours) {
  srat_neighbour <- FindNeighbors(srat, return.neighbor = TRUE)
  nn_obj <- srat_neighbour@neighbors$RNA.nn
  nn_cells <- nn_obj@nn.idx
  nn_dist <- nn_obj@nn.dist
  nn_id <- nn_obj@cell.names
  write.table(nn_cells, file=paste(output_path, "nn_cells.txt", sep=""), row.names=TRUE, col.names=FALSE, sep="\t", quote = FALSE)
  write.table(nn_dist, file=paste(output_path, "nn_dist.txt", sep=""), row.names=TRUE, col.names=FALSE, sep="\t", quote = FALSE)
  write.table(nn_id, file=paste(output_path, "nn_id.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
}
srat_graphs <- FindNeighbors(srat)
if (!do_neighbours) {
  nn_id <- srat_graphs@graphs$RNA_nn@Dimnames[[1]]
  write.table(nn_id, file=paste(output_path, "nn_id.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
}
nn_matrix <- as(as.matrix(srat_graphs@graphs$RNA_nn), "dgCMatrix")
snn_matrix <- as(as.matrix(srat_graphs@graphs$RNA_snn), "dgCMatrix")
writeMM(nn_matrix, paste(output_path, "nn.mtx", sep=""))
writeMM(snn_matrix, paste(output_path,"snn.mtx", sep=""))
# Algorithms in Seurat clustering:
# 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm
res_vect = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5) 
for ( alg in 1:4) {
  for ( res in res_vect) {
    outfile = paste("srat_a", as.character(alg), "_r", gsub("\\.", "", sprintf("%.1f", res)), sep="")
    srat_clust <- FindClusters(srat_graphs, resolution = res, algorithm = alg)
    write.table(srat_clust@active.ident, file=paste(output_path, outfile, sep=""), row.names=TRUE, col.names=FALSE, sep="\t", quote = FALSE)
  }
}
