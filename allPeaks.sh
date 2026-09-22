#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --partition=standard
#SBATCH --mem 7G
#SBATCH --account=scjp1
#SBATCH --output=%u-%x-%j.log
#SBATCH --error=%u-%x-%j.err
#SBATCH --mail-user=kgmoo@umich.edu
#SBATCH --mail-type=END,FAIL

Rscript /nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/sceAllPeaks.R
