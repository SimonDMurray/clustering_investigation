library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
theargs <- R.utils::commandArgs(asValues=TRUE)
input_seur <- theargs$seur
input_rcl <- theargs$rcl
input_bar <- theargs$barcodes
input_clus <- theargs$level
input_single <- theargs$single
output_path <- theargs$output
#Reading in Seurat object and copying data matrix to count matrix
h5_srat <- LoadH5Seurat(input_seur)
h5_srat@assays$RNA@counts <- h5_srat@assays$RNA@data
srat <- FindVariableFeatures(h5_srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
srat_graphs<- FindNeighbors(srat)
#Reading in rcl table
rcl_table <- read.table(input_rcl, sep = "\t", header = TRUE)
#Set nesting column variable type to factor
rcl_table$nesting <- as.factor(rcl_table$nesting)
#Reading in barcodes table
barcodes_table <- read.table(input_bar, header = TRUE)
colnames(barcodes_table) <- c("node", "barcode")
#Adding nesting column to metadata
srat_graphs@meta.data$nesting <- NA
#Selecting specifc cluster levels
cluster_levels = c()
list_levels <- str_split(input_clus, ",")
for (i in list_levels[[1]]) {
  cluster_levels <- append(cluster_levels, as.numeric(i))
}
level_vec = c()
#Selecting multiple levels to cluster 
for (row_num in 1:nrow(rcl_table)) {
  if (rcl_table[row_num, 1] %in% cluster_levels) {
    str_nodes <- rcl_table[row_num, 5]
    vec_nodes <- as.vector(strsplit(str_nodes, split = " "))[[1]]
    barcodes <- barcodes_table$barcode[barcodes_table$node %in% vec_nodes]
    srat_graphs@meta.data[barcodes, ]$nesting <- as.character(rcl_table[row_num, 4])
    add_nodes <- vec_nodes[!(vec_nodes %in% level_vec)]
    level_vec <- append(level_vec, add_nodes)
    level_bar <- barcodes_table$barcode[barcodes_table$node %in% level_vec]
  }
}
#Counting different levels of nesting in metadata
dplyr::count(srat_graphs@meta.data, nesting)
#Plotting clusterings with default colours
#Add "label = TRUE" to DimPlot to get nestings listed on plot
dp_default <- DimPlot(srat_graphs, group.by = "nesting") + NoLegend() + ggtitle("Clusters")
#Dimplot does have DiscretePalette options using pals: https://kwstat.github.io/pals/ but not enough colours
#Looked here for colour codes: https://www.designwizard.com/blog/design-trends/colour-combination
#Generate custom colouring
colours <- c("#195190FF", "#F4DF4EFF", "#FC766AFF", "#5B84B1FF", "#5F4B8BFF", "#E69A8DFF", 
             "#42EADDFF", "#CDB599FF", "#00A4CCFF", "#F95700FF", "#00203FFF", "#ADEFD1FF", 
             "#D6ED17FF", "#ED2B33FF", "#2C5F2DFF", "#EEA47FFF", "#9CC3D5FF", "#D198C5FF", 
             "#E0C568FF", "#CBCE91FF", "#EA738DFF", "#B1624EFF", "#5CC8D7FF", "#E3CD81FF",
             "#F2AA4CFF", "#A07855FF", "#D4B996FF", "#195190FF", "#603F83FF", "#CE4A7EFF", 
             "#2BAE66FF", "#FAD0C9FF", "#E94B3CFF", "#DAA03DFF", "#990011FF", "#435E55FF", 
             "#D64161FF", "#4B878BFF", "#76528BFF", "#FAEBEFFF", "#F93822FF", "#FDD20EFF", 
             "#755139FF", "#990011FF", "#435E55FF", "#CBCE91FF", "#333D79FF", "#F2EDD7FF",
             "#006B38FF", "#F95700FF", "#FFD662FF", "#00539CFF", "#D7C49EFF", "#343148FF", 
             "#FFA177FF", "#F5C7B8FF", "#DF6589FF", "#3C1053FF", "#FFE77AFF", "#2C5F2DFF",
             "#DD4132FF", "#9E1030FF", "#FCF951FF", "#422057FF", "#D01C1FFF", "#00B1D2FF")
#Plotting custom colouring
dp_custom <- DimPlot(srat_graphs, group.by = "nesting", cols = colours) + NoLegend() + ggtitle("Custom Clusters")
#Plotting specific clustering level
dp_level <- DimPlot(srat_graphs, group.by = "nesting", cells.highlight = level_bar, cols.highlight = "#004225", sizes.highlight = 0.1, label = TRUE) + NoLegend() + ggtitle(paste("Cluster level:", paste(cluster_levels, sep = ",")))
#Generates dictionary of all level sets of nodes where key is nesting (as a character) and value is list of nodes (as comma separated character)
single_dict <- vector(mode="list")
for (row_num in 1:nrow(rcl_table)) {
  single_2_vec = c()
  str_nodes <- rcl_table[row_num, 5]
  vec_nodes <- as.vector(strsplit(str_nodes, split = " "))[[1]]
  add_nodes <- vec_nodes[!(vec_nodes %in% single_2_vec)]
  single_2_vec <- append(single_2_vec, add_nodes)
  single_2_bar <- barcodes_table$barcode[barcodes_table$node %in% single_2_vec]
  bar_str <- paste(single_2_bar, collapse = ",")
  single_dict[as.character(rcl_table[row_num, 4])] <- bar_str
}
#Identify different nestings in dictionary
names(single_dict)
#Plotting single clustering at a specific level
bar_highlight <- strsplit(as.character(single_dict[input_single]), split = ",")[[1]]
dp_single <- DimPlot(srat_graphs, group.by = "nesting", cells.highlight = bar_highlight, cols.highlight = "#004225", sizes.highlight = 0.1) + NoLegend() + ggtitle(input_single)
#Export metadata table
metatable <- srat_graphs@meta.data
write.table(metatable, paste(output_path, "metadata.txt", sep = ""), sep = "\t", quote = FALSE)
#Identifying cell type and subtypes
factor_table_3_2 <- setDT(metatable)[, list(Subfactors = paste(unique(celltype.l3), collapse = " ")), by = celltype.l2]
colnames(factor_table_3_2) <- c("Type", "Subtype")
factor_table_2_1 <- setDT(metatable)[, list(Subfactors = paste(unique(celltype.l2), collapse = " ")), by = celltype.l1]
colnames(factor_table_2_1) <- c("Type", "Subtype")
#Exporting subtype tables
write.table(factor_table_3_2, paste(output_path, "factor_table_3_2.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(factor_table_2_1, paste(output_path, "factor_table_2_1.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
output_r = paste(output_path, "azimuth_vitessce_clustering.Rdata", sep = "")
save.image(output_r)