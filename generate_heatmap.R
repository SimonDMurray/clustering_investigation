library(gplots)
theargs <- R.utils::commandArgs(asValues=TRUE)
do_clus <- !is.null(theargs$rds)
output_path <- theargs$output
input_path <- theargs$input
mat <- read.table(input_path, header= TRUE)
names_of_rows <- mat$X0
mat <- mat[, -1]
rownames(mat) <- names_of_rows
mat <- as.matrix(mat)
transposed_mat <- t(mat)
# PLEASE NOTE: to execute on a cluster, the changing of plotting scale, axes margins etc has been removed. 
# Load into Rstudio and add the following settings to tinker with plot (feel free to change these settings if you prefer):
#
# png("output_name.png", width=1000, height=1000)
# heatmap.2( ... , cexRow=0.5, cexCol = 0.5, lhei=c(5, 30), lwid = c(5, 30)), margins = c(5, 5))
# dev.off()
#
if (do_clus) {
  hm <- heatmap.2(transposed_mat, trace = "none")
} else {
  hm <- heatmap.2(mat, trace = "none", Rowv = NULL, Colv = NULL)
}
output_r = paste(output_path, "heatmap.Rdata", sep = "")
save.image(output_r)
