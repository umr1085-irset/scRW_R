#!/bin/bash

#SBATCH --job-name="dag"
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

. /local/env/envconda.sh
conda activate renv
rm -r .snakemake/locks/
mkdir -p OUTPUT/DAG
snakemake --dag | dot -Tsvg > OUTPUT/DAG/dag.svg