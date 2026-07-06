#!/bin/bash

# Taking chopped up SingleCellAssay files and turning them into sex specific MAST runs looking for DEGs
# Max out at 180G for standard partition

#SBATCH -J MAST_SM_RNA
#SBATCH --account=scjp1
#SBATCH -t 24:00:00
#SBATCH --mem 256G
#SBATCH --output=%u-%j.log
#SBATCH --error=%u-%j.err
#SBATCH --mail-user=kgmoo@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=largemem

Rscript /nfs/turbo/umms-scjp/kgmoo/analysis/20240814_kgmoo_RNAge/RNA_MAST_T1.R
