#!/bin/bash

#SBATCH --job-name="snek"
#SBATCH --output=output_snek.out
#SBATCH --mem=400G
#SBATCH --cpus-per-task=50
#SBATCH --partition=sihp # cluster partition to use
#SBATCH --mail-user=<user-email> # user email
#SBATCH --mail-type=ALL # receive emails for all updates

. /local/env/envconda.sh
conda activate renv # activate R environment
rm -r .snakemake/locks/ # remove potential locks
snakemake --dag | dot -Tsvg > OUTPUT/DAG/dag.svg # generate workflow DAG
snakemake --resources load=5 -p -j 48 # launch
