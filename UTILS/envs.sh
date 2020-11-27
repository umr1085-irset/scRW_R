#!/bin/bash

#conda env remove --name renv
conda create -n renv python=3.7 r-base=3.6 --yes # create conda env with R
conda activate renv # activate env
#conda install -c conda-forge r=3.4.1 # install desired version of R
conda install -c bioconda bioconductor-scater --yes
conda install -c bioconda bioconductor-scran --yes
conda install -c conda-forge -c bioconda snakemake --yes
conda install -c bioconda bioconductor-dropletutils --yes
conda install -c bioconda bioconductor-hdf5array  --yes
# add rm HDF5Array
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("HDF5Array")
#BiocManager::install("scran")
conda install -c conda-forge r-robustbase --yes
conda install -c r r-kernsmooth --yes
conda install -c anaconda pytables --yes
# copy libRblas.so to .conda/envs/env/lib/R/lib/
#conda install -c conda-forge r-devtools
#devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
#conda install -c conda-forge r-remotes --yes
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#conda install -c bioconda r-seurat --yes 
conda install -c conda-forge umap-learn --yes
#remotes::install_github('ChristophH/sctransform')
#BiocManager::install("lfa")
conda install -c bioconda r-jackstraw --yes
conda install -c anaconda jinja2 --yes
conda install -c plotly plotly --yes
conda install -c r r-rjson --yes
conda install -c bokeh bokeh --yes