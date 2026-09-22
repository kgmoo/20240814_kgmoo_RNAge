##################
### Keagan Moo ###
### 2025/10/01 ###
### dbGap Expo ###
##################

# This script loads in and filters the ATAC and RNA total cross cell type counts matricies for the FUSION
# Skeletal Muscle data for samples that consented to data sharing. Subsequently the sample IDs are anonymized
# according to the common reference at:
# /nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/dbgap/molecular/sn_muscle_nucleus_cluster_sample_mapping.tsv

suppressPackageStartupMessages({library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(ggplot2)
library(glue)
})

CellTypes <- c("Type_1", "Type_2a", "Type_2x", "Endothelial", "Mesenchymal_Stem_Cell", "Macrophage", "Satellite_Cell", 
               "Neuromuscular_junction", "Neuronal", "Smooth_Muscle", "Muscle_Fiber_Mixed", "T_cell", "Adipocyte")
reFrame <- read.table("/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/dbgap/molecular/sn_muscle_nucleus_cluster_sample_mapping.tsv", sep = "\t", header = TRUE)

for (cDex in CellTypes) {
    
    print(glue("Loading and filtering {cDex} data..."))
    # Load ATAC
    small <- readRDS(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/{cDex}_FUSIONsnATAC_sce.rds"))
    smallF <- small[,small$sex == "F"]
    smallM <- small[,small$sex == "M"]
    keepRows <- rowSums(counts(small) >= 1) >= ncol(counts(small))*0.1
    small <- small[keepRows,]
    keepRowsF <- rowSums(counts(smallF) >= 1) >= ncol(counts(smallF))*0.1
    smallF <- smallF[keepRowsF,]
    keepRowsM <- rowSums(counts(smallM) >= 1) >= ncol(counts(smallM))*0.1
    smallM <- smallM[keepRowsM,]
    
    tempNucleiATAC <- dim(small)[2]
    
    # Load RNA
    sca <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nTotal_CDR_sca.rds"))
    scaF <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nFemale_CDR_sca.rds"))
    scaM <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nMale_CDR_sca.rds"))
    tempNucleiRNA <- dim(sca)[2]
    
    tDex <- "CleanUrFuckingData"
    
    if(cDex == "Mesenchymal_Stem_Cell") {
            tDex <- "FAP"
    } else if (cDex == "Neuromuscular_junction") {
            tDex <- "NMJ"
    } else {
            tDex <- cDex
    }
    
    tempreFrame <- reFrame[reFrame$coarse_cluster_name == glue("{tDex}"),]
    
    print(glue("ATAC Starting: {tempNucleiATAC} nuclei"))
    smaller <- small[,small$index %in% tempreFrame[tempreFrame$modality == "atac",]$index]
    tempNucleiATAC <- dim(smaller)[2]
    print(glue("ATAC Complete: {tempNucleiATAC} nuclei remain"))
    smallerF <- smallF[,smallF$index %in% tempreFrame[tempreFrame$modality == "atac",]$index]
    smallerM <- smallM[,smallM$index %in% tempreFrame[tempreFrame$modality == "atac",]$index]

    
    print(glue("RNA Starting: {tempNucleiRNA} nuclei"))
    scaler <- sca[,sca$id %in% tempreFrame[tempreFrame$modality == "rna",]$index]
    tempNucleiRNA <- dim(scaler)[2]
    print(glue("RNA Complete: {tempNucleiRNA} nuclei"))
    scalerF <- scaF[,scaF$id %in% tempreFrame[tempreFrame$modality == "rna",]$index]
    scalerM <- scaM[,scaM$id %in% tempreFrame[tempreFrame$modality == "rna",]$index]
    
    
    tempATACFrame <- data.frame(as.matrix(counts(smaller)))
    tempATACFrameF <- data.frame(as.matrix(counts(smallerF)))
    tempATACFrameM <- data.frame(as.matrix(counts(smallerM)))

    tempRNAFrame <- data.frame(as.matrix(counts(scaler)))
    tempRNAFrameF <- data.frame(as.matrix(counts(scalerF)))
    tempRNAFrameM <- data.frame(as.matrix(counts(scalerM)))

    write.table(tempATACFrame, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/snCountsMatrix_1ct10pct/FUSION_sampleTotal_tissueSM_nuclei{tDex}_atac1ct10pctN.tsv"), sep = "\t")
    write.table(tempATACFrameF, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/snCountsMatrix_1ct10pct/FUSION_sampleFemale_tissueSM_nuclei{tDex}_atac1ct10pctN.tsv"), sep = "\t")
    write.table(tempATACFrameM, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/snCountsMatrix_1ct10pct/FUSION_sampleMale_tissueSM_nuclei{tDex}_atac1ct10pctN.tsv"), sep = "\t")
    
    write.table(tempRNAFrame, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/snCountsMatrix_1ct10pct/FUSION_sampleTotal_tissueSM_nuclei{tDex}_rna1ct10pctN.tsv"), sep = "\t")
    write.table(tempRNAFrameF, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/snCountsMatrix_1ct10pct/FUSION_sampleFemale_tissueSM_nuclei{tDex}_rna1ct10pctN.tsv"), sep = "\t")
    write.table(tempRNAFrameM, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/snCountsMatrix_1ct10pct/FUSION_sampleMale_tissueSM_nuclei{tDex}_rna1ct10pctN.tsv"), sep = "\t")
}