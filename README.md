# scRW_R
Single Cell RNA-seq Workflow in R

## Pipeline presentation

## Installation
On a local computer or remote platform with `Git` installed, run:
```
git clone https://github.com/umr1085-irset/scRW_R.git
```

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
  
Additional sections in the config file should not be altered.

## Run pipeline

For a dry-run (simulated run with no computation performed), run the `snakemake` command with the `-n` option. The `-j` option is used to set the maximum number of cores to use. If not specified, all avaible cores are used.
```
snakemake -n -j 4
```

If the dry-run looks good, run the command without the `-n` option:
```
snakemake -j 4
```

The following example shows how to launch the workflow on a cluster using `slurm`. Users can use the `sbatch` command on a bash script `run_snakemake.sh` that looks like the following:
```
#!/bin/bash

#SBATCH --job-name="scrw" # job name
#SBATCH --output=output_scrw.out # job output file
#SBATCH --mem=200G # job memory
#SBATCH --cpus-per-task=48 # job cores
#SBATCH --partition=sihp # cluster partition to use
#SBATCH --mail-user=<e-mail> # user email
#SBATCH --mail-type=ALL # receive emails for all updates

. /local/env/envconda.sh # load Conda
conda activate renv # load appropriate R environment
rm -r .snakemake/locks/ # remove potential locks on the output folder
snakemake --dag | dot -Tsvg > OUTPUT/DAG/dag.svg # generate directed acyclic graph of workflow
snakemake --resources load=3 -p -j 48 # run workflow
```

To launch the workflow, use the following:
```
sbatch run_snakemake.sh
```
