##################
### Keagan Moo ###
### 2026/02/20 ###
### Motif Stat ###
##################

# Runs TF motif enrichment for age-associated peaks in defined cell type(s). 

suppressPackageStartupMessages({library(tidyverse)
library(matrixStats)
library(dplyr)
library(glue)
library(data.table)
library(tidyr)
library(hash)
})

Types <- c("Type_2x")

for (cDex in Types) {
    print(cDex)
    print("Starting ...")
    print("Loading metadata")
    metadata.factors <- readRDS("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/KMScripts/metadata.factors")
    peaks_dir <- file.path("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/",
                           "arushivData20220415/fusion_150_ccc")
    print("Loading background peaks")
    peaks <- file.path(peaks_dir, glue("atac.{cDex}.ccc.tsv")) %>%
      fread(data.table = F) %>%
      tibble::column_to_rownames("feature") %>%
      data.matrix()
    
    peaks <- peaks[, colnames(peaks) != "V1"]
    
    # Put everything in the same order and remove samples with missing age
    peaks <- peaks[,colnames(peaks) %in% rownames(metadata.factors)]
    metadata.factors <- metadata.factors[metadata.factors$sample %in% colnames(peaks),]
    peaks <- peaks[,rownames(metadata.factors)]
    peaksBed <- rownames(peaks)
    peaksBedDF <- t(data.frame(strsplit(peaksBed, split = ":")))
    rownames(peaksBedDF) <- NULL
    colnames(peaksBedDF) <- NULL
    
    print("Loading TF Set")
    pListTFs <- data.frame(data = list.files(path = "/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/chromatinTFAll/ref/"))
    pListTFs$data <- substr(pListTFs$data,1,nchar(pListTFs$data)-4)
    pListTFs$data <- gsub("::", "_and_", pListTFs$data)
    pListTFs$data <- gsub("-", "_dash_", pListTFs$data)
    pListTFs <- pListTFs$data
    combinedModels <- data.frame(matrix(0, 2, 6))
    colnames(combinedModels) <- c("F GLM", "M GLM", 
                    "F GLM UP", "M GLM UP", "F GLM DOWN", "M GLM DOWN")
    rownames(combinedModels) <- c("Test Coef", "Test PNom")
    
    peaksEns <- read.table(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/enrich/{cDex}_bgoEns.tsv"), sep = "\t")[,c("V1", "V2", "V3", "V10")]
    colnames(peaksEns) <- c("peak_chromosome", "peak_startPosition", "peak_endPosition", "Distance")
    peaksEns$Distance <- abs(peaksEns$Distance)
    peaksEns <- peaksEns %>% distinct(peak_chromosome, peak_startPosition, peak_endPosition, Distance, .keep_all=TRUE)

    # Load GLMNB defined age-associated peaks
    GLMNBSelectedF <- read.table(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/results/lmm/Agg/F_pTSS_Agg_{cDex}.tsv"), sep = "\t", header = TRUE)
    GLMNBSelectedM <- read.table(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/results/lmm/Agg/M_pTSS_Agg_{cDex}.tsv"), sep = "\t", header = TRUE)

    fComp <- GLMNBSelectedF %>% dplyr::select("estimate", "p_FDR") %>% dplyr::filter(p_FDR < 0.05)
    fCompUP <- fComp %>% dplyr::filter (estimate > 0)
    fCompDOWN <- fComp %>% dplyr::filter (estimate < 0)
    mComp <- GLMNBSelectedM %>% dplyr::select("estimate", "p_FDR") %>% dplyr::filter(p_FDR < 0.05)
    mCompUP <- mComp %>% dplyr::filter (estimate > 0)
    mCompDOWN <- mComp %>% dplyr::filter (estimate < 0)

    fSelectGLM <- data.frame(Select = rownames(fComp))
    fSelectGLMUP <- data.frame(Select = rownames(fCompUP))
    fSelectGLMDOWN <- data.frame(Select = rownames(fCompDOWN))
    mSelectGLM <- data.frame(Select = rownames(mComp))
    mSelectGLMUP <- data.frame(Select = rownames(mCompUP))
    mSelectGLMDOWN <- data.frame(Select = rownames(mCompDOWN))
    fpeaksSelectedGLM <- unique(fSelectGLM$Select)
    fpeaksSelectedGLMUP <- unique(fSelectGLMUP$Select)
    fpeaksSelectedGLMDOWN <- unique(fSelectGLMDOWN$Select)
    mpeaksSelectedGLM <- unique(mSelectGLM$Select)
    mpeaksSelectedGLMUP <- unique(mSelectGLMUP$Select)
    mpeaksSelectedGLMDOWN <- unique(mSelectGLMDOWN$Select)

    # Load in peaks that have a p-value for age-association greater than 0.5, we will use these as the background or "non age-associated set"
    AgePeaks <- readRDS(file = glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/KMScripts_Update/Res/ChromEnrich/AgePeaks_{cDex}_GLMMNB_05.rds"))

    nucData <- read.table(file = glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/enrich/{cDex}_Nuc.tsv"), sep = "\t")
    colnames(nucData) <- c("chr", "start", "end", "ATp", "GCp", "Act", "Cct", "Gct", "Tct", "Nct", "Otherct", "Total_Length")
    nucData$peaKey <- paste(nucData$chr, nucData$start, 
                            nucData$end, sep=":")
    GCData <- nucData %>% dplyr::select(c("peaKey", "GCp"))
    
    print("Beginning model loop")
    removedVec <- c()
    for (tempTF in pListTFs) {
        print(tempTF)
        colonTF <- gsub("_and_", "::", tempTF)
        colonTF <- gsub("_dash_", "-", colonTF)
        if (file.size(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/chromatinTFAll/intersect/{cDex}_bgoChrom_{colonTF}.tsv")) == 0L) {
            print("Empty File")
            peaksDF <- read.table(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/chromatinTFAll/removedForOverlap/{cDex}_bgoChrom_ARNT_3.tsv"), sep = "\t")
        } else {
            peaksDF <- read.table(glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20240124_kgmoo_WolverineScripts/Files/ARIMAke/chromatinTFAll/intersect/{cDex}_bgoChrom_{colonTF}.tsv"), sep = "\t")
        }
        if (ncol(peaksDF) == 9) {
        colnames(peaksDF) <- c("peak_chromosome", "peak_startPosition", "peak_endPosition", "tf_chromosome", "tf_startPosition", "tf_endPosition", "TF", "TF_affinity", "TF_strand")
        peaksDF$TF <- gsub("::", "_and_", peaksDF$TF)
        peaksDF$TF <- gsub("-", "_dash_", peaksDF$TF)
        peaksDF <- left_join(peaksDF, peaksEns, by=c("peak_chromosome", "peak_startPosition", "peak_endPosition"))
        libsizes <- colSums(peaks)
        peaksizes <- data.frame(peaKey=rownames(peaks), median=rowMedians(peaks))
        
        if (all(names(libsizes) == metadata.factors$sample)) {
            message("Adding library sizes to meta data")
            metadata.factors$libsize <- log2(unname(libsizes))
        } else {
            stop("Something wrong (2)")
        }
    
        peaKey <- paste(peaksDF$peak_chromosome, peaksDF$peak_startPosition, 
                    peaksDF$peak_endPosition, sep=":")
        inputDF <- data.frame(peaKey)
    
        inputDF$fSelectedGLM <- as.integer(inputDF$peaKey %in% fpeaksSelectedGLM)
        inputDF$mSelectedGLM <- as.integer(inputDF$peaKey %in% mpeaksSelectedGLM)
        inputDF$fSelectedGLMUP <- as.integer(inputDF$peaKey %in% fpeaksSelectedGLMUP)
        inputDF$mSelectedGLMUP <- as.integer(inputDF$peaKey %in% mpeaksSelectedGLMUP)
        inputDF$fSelectedGLMDOWN <- as.integer(inputDF$peaKey %in% fpeaksSelectedGLMDOWN)
        inputDF$mSelectedGLMDOWN <- as.integer(inputDF$peaKey %in% mpeaksSelectedGLMDOWN)
        inputDF$Distal5kb <- as.integer(peaksDF$Distance > 5000)
        inputDF <- merge(inputDF, peaksizes, by = "peaKey")
        inputDF <- merge(inputDF, GCData, by = "peaKey")
    
        inputDF <- inputDF[!duplicated(inputDF), ]
        duPeaks <- inputDF[duplicated(inputDF$peaKey),]$peaKey
        
        stateTerms <- peaksDF[peaksDF$TF != ".",c(7,1,2,3)]
    
        peaKollapse <- paste(stateTerms$peak_chromosome, stateTerms$peak_startPosition, 
                            stateTerms$peak_endPosition, sep=":")
        stateTerms$peakCoord <- peaKollapse
        statePeaks <- stateTerms[,c(1,5)]
        stateTerms <- statePeaks %>% group_by(TF) %>% summarise(peakCoord = list(peakCoord))
    
        FauxRaux <- c(0,1)
        setList <- c("F GLM", "M GLM", 
                    "F GLM UP", "M GLM UP", "F GLM DOWN", "M GLM DOWN")
        models <- data.frame(replicate(length(setList), FauxRaux, simplify=FALSE))
        colnames(models) <- setList
        pListStateC <- paste0(tempTF, " Coef")
        pListStateP <- paste0(tempTF, " PNom")
        # Interleaves the strings and then sets them to row names
        rownames(models) <- c(rbind(pListStateC,pListStateP))
    
        tryModel <- function(FUNformula, FUNinput, ...) {
           tempModel <- tryCatch(
             expr = {
               tempModel <- glm(formula = FUNformula, family="binomial", data=FUNinput)
             }, 
             error = function(e){
               print(glue("{tempName} suffered a fatal error"))
               tempModel <- -1
             }, 
             warning = function(w){
               print(glue("{tempName} suffered a warning"))
               tempModel <- -1
             }, 
             finally = {
               #pass
             }) 
           return(tempModel)
        }
    
        for (sDex in 1:nrow(stateTerms)){
      
          tempName <- stateTerms[sDex,1]
          print(tempName)
          inputDF[glue("{tempName}")] <- as.integer(inputDF$peaKey %in% noquote(stateTerms$peakCoord[[sDex]]))
     
        }
    
        setSet <- c("fSelectedGLM", "mSelectedGLM", 
                "fSelectedGLMUP", "mSelectedGLMUP", 
                "fSelectedGLMDOWN", "mSelectedGLMDOWN")
    
        peakSetSet <- list(fpeaksSelectedGLM, mpeaksSelectedGLM,
                           fpeaksSelectedGLMUP, mpeaksSelectedGLMUP,
                           fpeaksSelectedGLMDOWN, mpeaksSelectedGLMDOWN)

        for (pDex in 1:length(setSet)) {
           tempPeak <- setSet[pDex]
           tempPeakSet <- peakSetSet[pDex]
           tempAgePeaks <- dplyr::setdiff(AgePeaks, tempPeakSet[[1]])
           tempInputDF <- inputDF[!(inputDF$peaKey %in% tempAgePeaks),]#
      
            tempForm <- as.formula(glue("{tempPeak}~median+Distal5kb+GCp+{tempTF}"))
            tempModel <- tryModel(FUNformula = tempForm, FUNinput = tempInputDF)
            if((typeof(tempModel) == "list")){
                models[2*sDex - 1, pDex] <- tempModel$coefficients[5]
                models[2*sDex, pDex] <- summary(tempModel)$coefficients[20]
            }
        }
    
        combinedModels <- rbind(combinedModels, models)
        } else {
        removedVec <- c(tempTF, removedVec)
        print("---666---")
        }
    }
    
    print("Loop Complete")
    if (length(removedVec) > 0){
        print("Rejcted TFs, no overlap")
        print(removedVec)
    }
    combinedModels <- combinedModels[-(1:2),]
    modelsDF <- data.frame(t(combinedModels))
        for (nDex in 1:ncol(modelsDF)) {modelsDF[,nDex] <- as.numeric(modelsDF[,nDex])}
    print("Saving")
    saveRDS(inputDF, file = glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/TFenrich/
                            {cDex}_TFEnrichmentInputDF_All3k_IndivMod_BGxGLM05_AGxGLM005.rds"))
    saveRDS(modelsDF, file = glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/TFenrich/
                            {cDex}_TFEnrichmentModels_All3k_IndivMod_BGxGLM05_AGxGLM005.rds"))
    
    print("Melting")
    modelsCoef <- modelsDF %>% dplyr::select(contains("Coef"))
    modelsNom <- modelsDF %>% dplyr::select(contains("PNom"))
    modelsMelt <- data.frame(TF = -1, peakSet = -1, Coef = -1, NomP = -1)
    
    for(mDex in 1:ncol(modelsCoef)) {
      
      metaFrame <- data.frame(TF = rep(gsub(".Coef", "", colnames(modelsCoef)[mDex]), nrow(modelsDF)), peakSet = rownames(modelsCoef))
      metaFrame$Coef <- modelsCoef[,mDex]
      metaFrame$NomP <- modelsNom[,mDex]
      modelsMelt <- rbind(modelsMelt, metaFrame)
      
    }
    
    modelsMelt <- modelsMelt[-1,]
    modelsMelt$BonP <- p.adjust(modelsMelt$NomP, method = "bonferroni", n=nrow(modelsMelt))
    modelsMelt$Type <- glue("{cDex}")
    write.csv(modelsMelt, file = glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/TFenrich/
                                        {cDex}_TFEnrichmentModels_All3k_IndivMod_BGxGLM05_AGxGLM005.csv"), row.names = FALSE, quote = FALSE)
}