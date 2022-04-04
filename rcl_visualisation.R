library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(msigdbr)
theargs <- R.utils::commandArgs(asValues=TRUE)
input_seurat <- theargs$seurat
input_rcl <- theargs$rcl
input_barcode <- theargs$barcode
output_path <- theargs$output
#Reading in Seurat object and copying data matrix to count matrix
h5_srat <- LoadH5Seurat(input_seurat)
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
barcodes_table <- read.table(input_barcode, header = TRUE)
colnames(barcodes_table) <- c("node", "barcode")
#Adding nesting column to metadata
srat_graphs@meta.data$nesting <- NA
#selecting specifc cluster levels
cluster_levels = c(1,2)
level_select = 1
level_vec = c()
#selecting multiple levels to cluster 
for (row_num in 1:nrow(rcl_table)) {
  if (rcl_table[row_num, 1] %in% cluster_levels) {
    str_nodes <- rcl_table[row_num, 5]
    vec_nodes <- as.vector(strsplit(str_nodes, split = " "))[[1]]
    barcodes <- barcodes_table$barcode[barcodes_table$node %in% vec_nodes]
    srat_graphs@meta.data[barcodes, ]$nesting <- as.character(rcl_table[row_num, 4])
    #selecting a single specific level to isolate
    if (rcl_table[row_num, 1] == level_select) {
      add_nodes <- vec_nodes[!(vec_nodes %in% level_vec)]
      level_vec <- append(level_vec, add_nodes)
      level_bar <- barcodes_table$barcode[barcodes_table$node %in% level_vec]
    }
  }
}
dplyr::count(srat_graphs@meta.data, nesting)
#Plotting clusterings with default colours
# add label = TRUE to DimPlot to get nestings listed on plot
default_clusters <- DimPlot(srat_graphs, group.by = "nesting") + NoLegend() + ggtitle("Clusters")
#Dimplot does have DiscretePalette options using pals: https://kwstat.github.io/pals/ but not enough colours
#looked here for colour codes: https://www.designwizard.com/blog/design-trends/colour-combination
#generated custom colouring
colours <- c("#195190FF", "#F4DF4EFF", "#FC766AFF", "#5B84B1FF", "#5F4B8BFF", "#E69A8DFF", "#42EADDFF", "#CDB599FF", 
             "#00A4CCFF", "#F95700FF", "#00203FFF", "#ADEFD1FF", "#D6ED17FF", "#ED2B33FF", "#2C5F2DFF", "#EEA47FFF",
             "#9CC3D5FF", "#D198C5FF", "#E0C568FF", "#CBCE91FF", "#EA738DFF", "#B1624EFF", "#5CC8D7FF", "#E3CD81FF",
             "#F2AA4CFF", "#A07855FF", "#D4B996FF", "#195190FF", "#603F83FF", "#CE4A7EFF", "#2BAE66FF", "#FAD0C9FF",
             "#E94B3CFF", "#DAA03DFF", "#990011FF", "#435E55FF", "#D64161FF", "#4B878BFF", "#76528BFF", "#FAEBEFFF",
             "#F93822FF", "#FDD20EFF", "#755139FF", "#990011FF",
             "#435E55FF", "#CBCE91FF", "#333D79FF", "#F2EDD7FF",
             "#006B38FF", "#F95700FF", "#FFD662FF", "#00539CFF",
             "#D7C49EFF", "#343148FF", "#FFA177FF", "#F5C7B8FF",
             "#DF6589FF", "#3C1053FF", "#FFE77AFF", "#2C5F2DFF",
             "#DD4132FF", "#9E1030FF", "#FCF951FF", "#422057FF",
             "#D01C1FFF", "#00B1D2FF")
custom_clusters <- DimPlot(srat_graphs, group.by = "nesting", cols = colours) + NoLegend() + ggtitle("Custom Clusters")
level_clusters <- DimPlot(srat_graphs, group.by = "nesting", cells.highlight = level_bar, cols.highlight = "#004225", sizes.highlight = 0.1, label = TRUE) + NoLegend() + ggtitle(paste("Cluster level:", level_select))
# generates dictionary of all level sets of nodes where key is nesting (as a character) and value is list of nodes (as comma 
# separated character)
rcl_dict <- vector(mode="list")
for (row_num in 1:nrow(rcl_table)) {
  single_2_vec = c()
  str_nodes <- rcl_table[row_num, 5]
  vec_nodes <- as.vector(strsplit(str_nodes, split = " "))[[1]]
  add_nodes <- vec_nodes[!(vec_nodes %in% single_2_vec)]
  single_2_vec <- append(single_2_vec, add_nodes)
  single_2_bar <- barcodes_table$barcode[barcodes_table$node %in% single_2_vec]
  bar_str <- paste(single_2_bar, collapse = ",")
  rcl_dict[as.character(rcl_table[row_num, 4])] <- bar_str
}
#creating dictionary of all cluster top marker gene dataframes
df_dict <- vector(mode="list")
for ( row_num in 1:nrow(rcl_table)) {
  bar_highlight <- unlist(strsplit(as.character(rcl_dict[rcl_table[row_num, 4]]), split = ","))
  subset_srat <- srat_graphs[, bar_highlight]
  nest_markers <- FindAllMarkers(subset_srat)
  df_dict[as.character(rcl_table[row_num, 4])] <- list(nest_markers)
}
#creating new table that will contain marker gene info for each nesting
gene_table <- rcl_table[,c(1:4)]
gene_table$genes <- NA
gene_table$enrichment <- NA
gene_table$go_description <- NA
#removing cluster info from nesting column leaving just cluster size in table
nesting_vec <- as.character(rcl_table$nesting)
for (index in 1:length(nesting_vec)) {
  nest_split <- unlist(strsplit(nesting_vec[index], "::"))
  for (clust_size in 1:length(nest_split)) {
    size <- gsub(".*_","",nest_split[clust_size])
    nest_split[clust_size] <- size
  }
  size_string <- paste(as.character(nest_split), collapse="_")
  nesting_vec[index] <- size_string
}
gene_table$nesting <- nesting_vec
# Adding top marker gene info to table
for (row_num in 1:nrow(rcl_table)) {  
  nesting <- as.character(rcl_table[row_num, 4])
  nest_markers <- as.data.frame(df_dict[nesting])
  if (length(nest_markers) == 0) {
    write.table(nest_markers, paste(output_path, "NO_MARKERS_", nesting, ".txt", sep = ""),  quote = FALSE, sep = "\t", row.names = FALSE)
  }
  if (length(nest_markers) != 0) {
    colnames(nest_markers) <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
    write.table(nest_markers, paste(output_path, nesting, ".txt", sep = ""),  quote = FALSE, sep = "\t", row.names = FALSE)
    log_p_vals <- -log10(nest_markers$p_val)
    nest_markers$rounded_p_val <- round(log_p_vals)
    nest_markers$rounded_p_val <- as.numeric(gsub("Inf", 323, nest_markers$rounded_p_val))
    nest_markers <- nest_markers %>% arrange(desc(rounded_p_val))
    top_genes <- head(rownames(nest_markers), n=10)
    top_pvals <- head(nest_markers$rounded_p_val, n=10)
    combined_vec <- c()
    for(index in 1:length(top_genes)){
      combined_char <- paste(top_genes[index], "(", top_pvals[index], ")", sep = "")
      combined_vec <- append(combined_vec, combined_char)
      combined_string <- paste(as.character(combined_vec), collapse=",")
    }
    gene_table[row_num, 5] <- combined_string
  }
}
#adding GO database info to table
for (row_num in 1:nrow(rcl_table)) {
  nesting <- as.character(rcl_table[row_num, 4])
  nest_markers <- as.data.frame(df_dict[nesting])
  if (length(nest_markers) == 0) {
    write.table(nest_markers, paste(output_path, "NO_MARKERS_", nesting, ".txt", sep = ""),  quote = FALSE, sep = "\t", row.names = FALSE)
  }
  if (length(nest_markers) != 0) {
    colnames(nest_markers) <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
    nest_genes <- nest_markers$gene
    nest_entrez <- mapIds(org.Hs.eg.db, nest_genes, 'ENTREZID', 'SYMBOL')
    no_nas <- as.vector(nest_entrez[!is.na(nest_entrez)])
    go_output <- enrichGO(no_nas, "org.Hs.eg.db") #, ont = "ALL")
    go_result <- go_output@result
    write.table(go_result, paste(output_path, nesting, ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
    go_result <- go_result %>% arrange(desc(p.adjust))
    top_go <- head(rownames(go_result), n=10)
    go_pvals <- round(head(go_result$p.adjust, n=10), digits = 4)
    go_desc <- head(go_result$Description, n=10)
    go_vec <- c()
    for(index in 1:length(top_go)){
      go_char <- paste(top_go[index], "(", go_pvals[index], ")", sep = "")
      go_vec <- append(go_vec, go_char)
      go_string <- paste(as.character(go_vec), collapse=",")
    }
    gene_table[row_num, 6] <- go_string
    gene_table[row_num, 7] <- paste(as.character(go_desc), collapse=", ")
  }
}
#exporting marker genes table
write.table(gene_table, paste(output_path, "markers.txt", sep = ""),  quote = FALSE, sep = "\t", row.names = FALSE)
#export metadata table
metatable <- srat_graphs@meta.data
write.table(metatable,  paste(output_path, "metadata.txt", sep = ""), sep = "\t", quote = FALSE)
#identifying cell type and subtypes
factor_table_3_2 <- setDT(metatable)[, list(Subfactors = paste(unique(celltype.l3), collapse = " ")), by = celltype.l2]
colnames(factor_table_3_2) <- c("Type", "Subtype")
factor_table_2_1 <- setDT(metatable)[, list(Subfactors = paste(unique(celltype.l2), collapse = " ")), by = celltype.l1]
colnames(factor_table_2_1) <- c("Type", "Subtype")
#exporting subtype tables
write.table(factor_table_3_2, paste(output_path, "factor_table_3_2.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(factor_table_2_1, paste(output_path, "factor_table_2_1.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
output_r = paste(output_path, "rcl_visualisation.Rdata", sep = "")
save.image(output_r)
