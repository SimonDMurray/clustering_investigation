library(Seurat)
theargs <- R.utils::commandArgs(asValues=TRUE)
output_path <- theargs$output
input_path <- theargs$input
do_rds <- !is.null(theargs$rds)
output_path <- theargs$output
input_path <- theargs$input
if (do_rds) {
  srat <- readRDS(input_path)
} else {
  matrix_data <- Read10X_h5(input_path, use.names = T)
  srat  <- CreateSeuratObject(counts = matrix_data)
}
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")
meta_data <- srat@meta.data
colnames(meta_data) <- c("orig.ident", "UMI_reads", "expressed_genes", "percent.mt", "percent.rb")
stats <- meta_data[c("UMI_reads", "expressed_genes", "percent.mt", "percent.rb")]
write.table(stats, file=paste(output_path, "bio_stats.txt", sep=""), row.names=TRUE, col.names = TRUE)
