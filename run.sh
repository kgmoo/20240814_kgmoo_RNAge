#!/bin/bash
#SBATCH -t 50:00:00
#SBATCH --mem 4000M
#SBATCH --account=scjp1
#SBATCH --output=%u-%x-%j.log
#SBATCH --error=%u-%x-%j.err
#SBATCH --mail-user=kgmoo@umich.edu
#SBATCH --mail-type=END,FAIL

module load singularity
workflow='glmmnb_FUSION'
nfdir='/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/scripts'
nextflow -C ${nfdir}/nextflow.config run ${nfdir}/main.nf -params-file config.yaml -resume -with-trace trace.txt 
