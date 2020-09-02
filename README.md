# scRW_R
Single Cell RNA-seq Workflow in R

## Pipeline presentation

## Config file
Before running the `snakemake` command to launch the analysis workflow, it is recommended to review the `config.yaml` file. It has multiple sections of interest:
* `DEFINE PATHS`:  
  This section is dedicated to paths and options relative to folders and files.
  * `OUTDIR`: output directory name
  * `AGGRMATRIX`: TESTDATA/small_aggr/filtered_feature_bc_matrix.h5 # aggregated H5 file
  * `AGGRFILE`: TESTDATA/small_aggr/aggregation.csv # aggregation file used with cellranger aggre
  * `SAMPLE_EXTRACTION_FROMCOL`: library_id # library_id or molecule_h5
  * `INDIVDIR`: TESTDATA/unaggr # individual samples directory


## Run pipeline
