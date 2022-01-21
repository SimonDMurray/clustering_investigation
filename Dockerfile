FROM rocker/r-ver:4.0.4
RUN apt-get update && apt-get install -y liblzma-dev libbz2-dev zlib1g libpng-dev libxml2-dev \
    gfortran-7 libglpk-dev libhdf5-dev libcurl4-openssl-dev img2pdf python3 python3-dev python3-pip

RUN pip3 install numpy==1.20 pandas umap-learn leidenalg igraph anndata regex
RUN Rscript -e "install.packages(c('BiocManager', 'devtools', 'R.utils'))"
RUN Rscript -e "BiocManager::install(c('Rhtslib', 'LoomExperiment', 'SingleCellExperiment'))"
RUN Rscript -e "devtools::install_github('cellgeni/sceasy')"
RUN Rscript -e "install.packages(c('hdf5r','dimRed','png','ggplot2','reticulate','plotly','Matrix','network','leiden', 'Seurat', 'gplots'))" 
