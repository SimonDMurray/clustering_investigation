library(Seurat)
theargs <- R.utils::commandArgs(asValues=TRUE)
do_rds <- !is.null(theargs$rds)
output_path <- theargs$output
input_path <- theargs$input
vol_path <- theargs$volatility
if (do_rds) {
  srat <- readRDS(input_path)
} else {
  adj.matrix <- Read10X(input_path)
  srat <- CreateSeuratObject(adj.matrix) 
}
DefaultAssay(srat) <- "ADT"
vol_file <- read.table(vol_path, sep="\t")
vol_df <- as.data.frame(vol_file$V2)
rownames(vol_df) <- vol_file$V1
srat[["volatility"]] <- vol_df
#srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat<- RunPCA(srat, features = VariableFeatures(object = srat))
srat <- FindNeighbors(srat)
srat <- FindClusters(srat)
srat <- RunUMAP(srat, dims = 1:50)
dp <- DimPlot(srat, label.size = 4,repel = T,label = T)
fp <- FeaturePlot(srat, features = "volatility", label = TRUE)
output_r = paste(output_path, "visualisation.Rdata", sep = "")
save.image(output_r)
