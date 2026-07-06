##################
### Keagan Moo ###
### 2024/10/09 ###
###  RNA MAST  ###
##################

# This script loads in the nuclei level data from the single nucleus skeletal muscle data from FUSION,
# and applies MAST to generate differential expression results tables
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("MAST"))
suppressPackageStartupMessages(library("lme4"))
suppressPackageStartupMessages(library("rsvd"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("glue"))

# Get cell type specific sca file
print("One")
clusterIO <- c("Type_1", "Type_2a", "Type_2x")

for (cDex in clusterIO){
sca <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nTotal_CDR_sca.rds"))
scaF <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nFemale_CDR_sca.rds"))
scaM <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nMale_CDR_sca.rds"))

# Load in complete metadata
print("Two")
metadataComp <- readRDS("/nfs/turbo/umms-scjp/kgmoo/data/skeletalmuscleRNA_metadata.rds") %>% dplyr::filter(coarse_cluster_name == cDex)
metadataCompF <- metadataComp[metadataComp$sex.x == "F",]
metadataCompM <- metadataComp[metadataComp$sex.x == "M",]

print("Three")

#Subset and sort colData as a more manageable dataframe, add ogtt metadata to each for use as a covariate
colOut <- colData(sca)
colOut <- data.frame(merge(colOut, metadataComp, by.x = "sample", by.y = "sample.x"))

colOut <- colOut %>% dplyr::arrange(id)
sca$ogtt <- colOut$ogtt

print("Four")

colOutF <- colData(scaF)
colOutF <- data.frame(merge(colOutF, metadataCompF, by.x = "sample", by.y = "sample.x"))
colOutF <- colOutF %>% dplyr::arrange(id)
scaF$ogtt <- colOutF$ogtt

print("Five")

colOutM <- colData(scaM)
colOutM <- data.frame(merge(colOutM, metadataCompM, by.x = "sample", by.y = "sample.x"))
print(dim(colOutM))
colOutM <- colOutM %>% dplyr::arrange(id)
scaM$ogtt <- colOutM$ogtt

print("Six")
# Hard set the datatypes for incoming covariates just in case
sca$sex <- as.factor(sca$sex)

sca$batch <- as.factor(sca$batch)
scaF$batch <- as.factor(scaF$batch)
scaM$batch <- as.factor(scaM$batch)

sca$ogtt <- as.factor(sca$ogtt)
scaF$ogtt <- as.factor(scaF$ogtt)
scaM$ogtt <- as.factor(scaM$ogtt)

sca$bmi <- as.numeric(sca$bmi)
scaF$bmi <- as.numeric(scaF$bmi)
scaM$bmi <- as.numeric(scaM$bmi)

sca$CDR <- as.numeric(sca$CDR)
scaF$CDR <- as.numeric(scaF$CDR)
scaM$CDR <- as.numeric(scaM$CDR)

sca$age <- as.numeric(sca$age)
scaF$age <- as.numeric(scaF$age)
scaM$age <- as.numeric(scaM$age)

print(head(colData(sca)))

# Run zlm models, strictConvergence and fitArgsD settings are from Maya Bose PhD which allow MAST to run in a resonable time frame
zlmCond_glm <- zlm(~ sex + batch + bmi + ogtt + CDR + (1|sample) + age, 
                   sca, 
                   method = 'glmer', 
                   ebayes = FALSE, 
                  exprs_value = 'logcounts',
                  strictConvergence = FALSE, 
                   fitArgsD = list(nAGQ = 0))
base::saveRDS(zlmCond_glm, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_zlm_CDR.rds"))

zlmCond_glmF <- zlm(~ batch + bmi + ogtt + CDR + (1|sample) + age, 
                   scaF, 
                   method = 'glmer', 
                   ebayes = FALSE, 
                  exprs_value = 'logcounts',
                  strictConvergence = FALSE, 
                   fitArgsD = list(nAGQ = 0))
summary(zlmCond_glmF)
base::saveRDS(zlmCond_glmF, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_Fzlm_CDR.rds"))

zlmCond_glmM <- zlm(~ batch + bmi + ogtt + CDR + (1|sample) + age, 
                   scaM, 
                   method = 'glmer', 
                   ebayes = FALSE, 
                  exprs_value = 'logcounts',
                  strictConvergence = FALSE, 
                   fitArgsD = list(nAGQ = 0))
summary(zlmCond_glmM)
base::saveRDS(zlmCond_glmM, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_Mzlm_CDR.rds"))

print("Seven")
zlmlr <- summary(zlmCond_glm, doLRT='age') 
zlmlrF <- summary(zlmCond_glmF, doLRT='age') 
zlmlrM <- summary(zlmCond_glmM, doLRT='age') 

zlmDT <- zlmlr$datatable
zlmDTF <- zlmlrF$datatable
zlmDTM <- zlmlrM$datatable

base::saveRDS(zlmDT, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_resDT_CDR.rds"))
base::saveRDS(zlmDTF, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_resDTF_CDR.rds"))
base::saveRDS(zlmDTM, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_resDTM_CDR.rds"))
print("Eight")
    
# Generate the hurdle model outputs and save them, technically this is easy to do outside of a SLURM job, but it is the ouput most people will want to see
# eg. one gene, one p-value, one coef

fcHurdle <- merge(zlmDT[contrast=='conditionStim' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      zlmDT[contrast=='conditionStim' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')

fcHurdleF <- merge(zlmDTF[contrast=='conditionStim' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      zlmDTF[contrast=='conditionStim' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdleF[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleFSig <- merge(fcHurdleF[fdr<.05 & abs(coef)>0], as.data.table(mcols(scaF)), by='primerid')

fcHurdleM <- merge(zlmDTM[contrast=='conditionStim' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      zlmDTM[contrast=='conditionStim' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdleM[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleMSig <- merge(fcHurdleM[fdr<.05 & abs(coef)>0], as.data.table(mcols(scaM)), by='primerid')

base::saveRDS(fcHurdle, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_hurDT_CDR.rds"))
base::saveRDS(fcHurdleF, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_hurDTF_CDR.rds"))
base::saveRDS(fcHurdleM, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_hurDTM_CDR.rds"))

}

print("one hundred and Done")