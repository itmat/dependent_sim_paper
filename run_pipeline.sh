#!/usr/bin/env sh
module load /project/itmatlab/sharedmodules/use.shared
module load python/3.10
module load R/4.3
source venv/bin/activate


### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --use-singularity --profile lsf -j 100 -c 100 "$@"
