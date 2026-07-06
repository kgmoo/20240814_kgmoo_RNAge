##################
### Keagan Moo ###
### 2024/10/09 ###
###  RNA MAST  ###
##################

# This script loads in the nuclei and sample level data from the single nucleus skeletal muscle data from FUSION,
# applies MILO to calculate relative nuclei cluster abundance and output both age-association 
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
cDex <- "Type_1"
#sca <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_10in10p_Drop5Samp_sca.rds"))
sca <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nTotal_CDR_sca.rds"))
scaF <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nFemale_CDR_sca.rds"))
scaM <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nMale_CDR_sca.rds"))

# Load in complete metadata, this step is out of sequence, so you'll need to generate this file after RNA_BREAK seperately
print("Two")
metadataComp <- readRDS("/nfs/turbo/umms-scjp/kgmoo/data/skeletalmuscleRNA_metadata.rds") %>% dplyr::filter(coarse_cluster_name == cDex)
metadataCompF <- metadataComp[metadataComp$sex.x == "F",]
metadataCompM <- metadataComp[metadataComp$sex.x == "M",]

# Load in somewhat redundant metadata to save on headache of bringing in per cluster live calculated metadata
#metadataFactors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/cluster_info_qc.tsv", sep = "\t", header = TRUE)
#batchdataFactors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/batch-sample-info-full.tsv", sep = "\t", header = TRUE)
#metadataFactors <- metadataFactors[(metadataFactors$modality == "rna") & (metadataFactors$demux_cat != "multiome"),]
#metadataFactors <- merge(metadataFactors, batchdataFactors, by = "SNG.1ST")

#metadataFactorsClust <- metadataFactors[(metadataFactors$coarse_cluster_name == cDex),]

# Live calculate the genes expressed within the cluster to each sample and nucleus
#scaLogCounts <- data.frame(logcounts(sca))
#colnames(scaLogCounts) <- colData(sca)$id
#scaNucleiGs <- data.frame(colSums(scaLogCounts > 0))
#rownames(scaNucleiGs) <- colData(sca)$id

#metadataJustAge <- metadataFactorsClust %>% select(index, sex, age, SNG.1ST, sample.x)
#metadataJustAge <- metadataJustAge[match(colnames(scaLogCounts), metadataJustAge$index),]
#metadataJustAge <- metadataJustAge[complete.cases(metadataJustAge),]
#rownames(metadataJustAge) <- metadataJustAge$index

#nucGs <- merge(scaNucleiGs, metadataJustAge, by = "row.names") 
#nucGs$Row.names <- NULL
#colnames(nucGs)[1] <- "Gs"

#print(head(nucGs))

#nucFGs <- nucGs[nucGs$sex == "F",] %>% select(index, Gs)
#nucMGs <- nucGs[nucGs$sex == "M",] %>% select(index, Gs)
#nucGs <- nucGs %>% select(index, Gs)

#print(dim(nucGs))
#print(dim(nucFGs))
#print(dim(nucMGs))

#scaLogCountS <- data.frame(t(scaLogCounts))
#metadataJustSamp <- metadataJustAge %>% select(sample.x)
#scaLogCountS <- merge(scaLogCountS, metadataJustSamp, by = "row.names")
#scaLogCountS$Row.names <- NULL

#print("Two point two")

#metadataJustSamp <- metadataJustAge %>% select(sample.x, sex, age) %>% unique()
#sampList <- unique(metadataJustSamp$sample.x)
#sampInit <- sampList[1]
#sampList <- sampList[-1]
##print(colnames(scaLogCountS))
#seedFrame <- (scaLogCountS[scaLogCountS$sample.x == sampInit,])
#seedFrame$sample.x <- NULL
##print(head(seedFrame))
#seedVec <- colSums(seedFrame)
#seedVal <- sum(seedVec != 0)
#valVec <- c(seedVal)

#print("Two point four")

#for (sampleDex in sampList) {

#    tempFrame <- (scaLogCountS[scaLogCountS$sample.x == sampleDex,])
#    tempFrame$sample.x <- NULL

#    tempVec <- colSums(tempFrame)
#    tempVal <- sum(tempVec != 0)
#    valVec <- c(valVec, tempVal)
    
#}

#print("Two point six")

#sampleGs <- data.frame(sample.x = metadataJustSamp$sample.x, age = metadataJustSamp$age, sex = metadataJustSamp$sex, Gs = valVec)
#sampleFGs <- sampleGs[sampleGs$sex == "F",] %>% select(sample.x, Gs)
#sampleMGs <- sampleGs[sampleGs$sex == "M",] %>% select(sample.x, Gs)
#sampleGs <- sampleGs %>% select(sample.x, Gs)

print("Two point eight")

#Subset and sort colData as a more manageable dataframe, add ogtt metadata to each for use as a covariate
colOut <- colData(sca)
colOut <- data.frame(merge(colOut, metadataComp, by.x = "sample", by.y = "sample.x"))
#colOut <- data.frame(merge(colOut, sampleGs, by.x = "sample", by.y = "sample.x"))
#colOut <- data.frame(merge(colOut, nucGs, by.x = "id", by.y = "index"))
colOut <- colOut %>% dplyr::arrange(id)
sca$ogtt <- colOut$ogtt
#sca$genesPS <- colOut$Gs
#sca$genesPN <- colOut$Gs

print("Two point nine")

colOutF <- colData(scaF)
colOutF <- data.frame(merge(colOutF, metadataCompF, by.x = "sample", by.y = "sample.x"))
#colOutF <- data.frame(merge(colOutF, sampleFGs, by.x = "sample", by.y = "sample.x"))
print(dim(colOutF))
#colOutF <- data.frame(merge(colOutF, nucFGs, by.x = "id", by.y = "index"))
colOutF <- colOutF %>% dplyr::arrange(id)
scaF$ogtt <- colOutF$ogtt
#scaF$genesPS <- colOutF$Gs
#scaF$genesPN <- colOutF$Gs

print("Two point ninenine")

colOutM <- colData(scaM)
colOutM <- data.frame(merge(colOutM, metadataCompM, by.x = "sample", by.y = "sample.x"))
#colOutM <- data.frame(merge(colOutM, sampleMGs, by.x = "sample", by.y = "sample.x"))
print(dim(colOutM))
#colOutM <- data.frame(merge(colOutM, nucMGs, by.x = "id", by.y = "index"))
colOutM <- colOutM %>% dplyr::arrange(id)
scaM$ogtt <- colOutM$ogtt
#scaM$genesPS <- colOutM$Gs
#scaM$genesPN <- colOutM$Gs

#colOld <- colData(sca)
#colOld <- data.frame(merge(colOld, metadataComp, by.x = "sample", by.y = "sample.x"))
#colOld <- colOld %>% dplyr::arrange(id)
#sca$ogtt <- colOld$ogtt

# Just a quick PCA, not that we don't trust prior handlers, but what if I want to check. Code adapted from MAST with RNASeq Vignette
#print("Check PCA")
#set.seed(123)
#plotPCA <- function(sca_obj){
#    projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
#    colnames(projection)=c("PC1","PC2","PC3","PC4")
#    pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
#    write.csv(pca , file = "/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/TestPCA.csv")
#    print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'PC4', 'sex', 'ogtt', 'age'),
#            mapping=aes(color=sex), upper=list(continuous='blank')))
#    invisible(pca)
#}

#pdf(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240814_kgmoo_RNAge/plots/PCA_SCACheck_{cDex}.pdf"), width=16, height=16)
#plotPCA(scaTPM)
#dev.off()

##### DO WE WANT TO EXPLORE ADAPTIVE THRESHOLDING? #####

#Subset the data by sex
#scaFOld <- subset(sca, sex == "F")
print("Three")
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

#sca$genesPS <- as.numeric(sca$genesPS)
#scaF$genesPS <- as.numeric(scaF$genesPS)
#scaM$genesPS <- as.numeric(scaM$genesPS)

#sca$genesPN <- as.numeric(sca$genesPN)
#scaF$genesPN <- as.numeric(scaF$genesPN)
#scaM$genesPN <- as.numeric(scaM$genesPN)

#sca$totalCounts <- as.numeric(sca$totalCounts)
#scaF$totalCounts <- as.numeric(scaF$totalCounts)
#scaM$totalCounts <- as.numeric(scaM$totalCounts)

sca$CDR <- as.numeric(sca$CDR)
scaF$CDR <- as.numeric(scaF$CDR)
scaM$CDR <- as.numeric(scaM$CDR)

sca$age <- as.numeric(sca$age)
scaF$age <- as.numeric(scaF$age)
scaM$age <- as.numeric(scaM$age)

print(head(colData(sca)))

# Run zlm models, strictConvergence and fitArgsD are from Maya Bose which allow MAST to run in a resonable time frame
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

print("Four")
zlmlr <- summary(zlmCond_glm, doLRT='age') 
zlmlrF <- summary(zlmCond_glmF, doLRT='age') 
#print(dimnames(zlmlrF))
zlmlrM <- summary(zlmCond_glmM, doLRT='age') 
#print(dimnames(zlmlrM))

zlmDT <- zlmlr$datatable
zlmDTF <- zlmlrF$datatable
zlmDTM <- zlmlrM$datatable

base::saveRDS(zlmDT, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_resDT_CDR.rds"))
base::saveRDS(zlmDTF, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_resDTF_CDR.rds"))
base::saveRDS(zlmDTM, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10n_resDTM_CDR.rds"))
print("Five")
# More code adapted from MAST vignette

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
print("one hundred and Done")
