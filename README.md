# scRW_R
Single Cell RNA-seq Workflow in R

## Pipeline presentation

## Config file
Before running the `snakemake` command to launch the analysis workflow, it is recommended to review the `config.yaml` file. It has multiple sections of interest:
* `DEFINE PATHS`:  
  This section is dedicated to paths and options relative to folders and files.
  * `OUTDIR`: output directory name
  * `AGGRMATRIX`: path to aggregated data H5 file
  * `AGGRFILE`:  aggregation file (.csv) used with cellranger aggr
  * `SAMPLE_EXTRACTION_FROMCOL`: `library_id` or `molecule_h5`. Users specify which column to use in the aggregation file to extract sample names
  * `INDIVDIR`: path to individual samples directory

* `WORKFLOW STEPS - set to True or False`:  
  Various steps are available in this section. Users can set them to `True` or `False` in order to generate their desired workflow. The latter will be generated automatically based on the binary values from this section.  
  
Additional sections below should not be altered.

## Run pipeline


