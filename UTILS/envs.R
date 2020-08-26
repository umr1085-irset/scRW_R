#add rm HDF5Array
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("HDF5Array")
BiocManager::install("scran")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
remotes::install_github('ChristophH/sctransform')
install.packages("rlang") # update rlang to 0.4.6
install.packages('vctrs') # update to 0.3.0
BiocManager::install("lfa")