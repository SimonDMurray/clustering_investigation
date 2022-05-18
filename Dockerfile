FROM rocker/r-ver:4.0.4
RUN apt-get update && apt-get install -y liblzma-dev libbz2-dev zlib1g libpng-dev libxml2-dev \
    gfortran-7 libglpk-dev libhdf5-dev libcurl4-openssl-dev img2pdf python3.8 python3-dev python3-pip

RUN pip3 install numpy==1.20.1 pandas umap-learn leidenalg igraph anndata regex
RUN Rscript -e "install.packages(c('BiocManager', 'devtools', 'R.utils', 'reticulate', 'processx'))"
RUN Rscript -e "BiocManager::install(c('Rhtslib', 'LoomExperiment', 'SingleCellExperiment', 'clusterProfiler', 'org.Hs.eg.db'))"
RUN Rscript -e "devtools::install_github(c('cellgeni/sceasy', 'satijalab/seurat-data'))"
RUN Rscript -e "install.packages(c('hdf5r','dimRed','png','ggplot2','reticulate','plotly','Matrix','network','leiden', 'Seurat', 'gplots', 'SoupX'))" 
RUN Rscript -e "install.packages('SeuratObject')"
RUN Rscript -e "devtools::install_github('mojaveazure/seurat-disk')"
