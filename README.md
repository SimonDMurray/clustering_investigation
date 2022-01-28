# clustering_investigation
Investigation into Clustering for comparison with MCL


biology-info.R - calculates percentage of ribosomal and mitochondrial DNA present

clustering.R - clusters Seurat data with Louivan, Multi-level Louvain, Smart Local Alogrithm and Leidan methods at various resolutiosn for clustering comparison. Generates matrix market format data of the nearest neighbour graphs.

combine_ID_weight.py - combines distance and cells neighbour files generated by Seurat for an alternative to martix market format.

convert_png.sh - converts Heatmaps saved as PNG files in R to PDF (R saves pdf files in a messy format).

generate_heatmap.R - generates heatmaps of clusterings and distance between each cluster in a pairwise format.

pca_statistics.R - normalises, finds variable features, scales and runs PCA on Seurat data before calcuating euclidean and dispersion values for each barcode.

produce_matrix.py - uses regex to obtain clustering and distandce data from file generated by clm dist and saves to tsv file.

run_r_scripts.sh - wrapper that allows user to execute any R script in this repo.

visualise_clustering.R - Seurat data that has been normalised, scaled and had UMAP ran on it is visualised to show clusterings and whether a node often changes cluster or not depending on method (mcxdump generates file used here)
