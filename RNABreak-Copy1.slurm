#!/bin/bash

# Get sample level gene count matricies

#SBATCH -J RNA_BREAK
#SBATCH --account=scjp1
#SBATCH -t 1:00:00
#SBATCH --mem 128G
#SBATCH --output=%u-%j.log
#SBATCH --error=%u-%j.err
#SBATCH --mail-user=kgmoo@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=largemem

Rscript /nfs/turbo/umms-scjp/kgmoo/analysis/20240814_kgmoo_RNAge/RNA_BREAK-Copy1.R
