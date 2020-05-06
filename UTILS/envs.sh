#!/bin/bash

#conda env remove --name renv
conda create -n renv python=3.7 r-base=3.6 --yes # create conda env with R
conda activate renv # activate env
#conda install -c conda-forge r=3.4.1 # install desired version of R
conda install -c bioconda bioconductor-scater --yes
conda install -c conda-forge -c bioconda snakemake --yes
conda install -c bioconda bioconductor-dropletutils --yes
conda install -c bioconda bioconductor-hdf5array  --yes
# add rm HDF5Array
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("HDF5Array")
#BiocManager::install("scran")
conda install -c conda-forge r-robustbase --yes
conda install -c r r-kernsmooth
conda install -c anaconda pytables --yes

#devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')