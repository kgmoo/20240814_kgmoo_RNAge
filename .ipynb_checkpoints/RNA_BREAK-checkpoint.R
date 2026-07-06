##################
### Keagan Moo ###
### 2024/08/15 ###
### RNA Break  ###
##################

# This script loads in the nuclei and sample level data from the single nucleus skeletal muscle data from FUSION post Arushi Varshney PhD QC,
# and breaks them up into filtered gene x count matricies on a per cell type basis
library("Seurat")
library("glue")
library("SingleCellExperiment")
library("scater")
library("scran")
library("tidyverse")
library("data.table")
library("MAST")

# Get list of rna count files as well as genes that occur in them from priot work done by Arushi Varshney PhD
print("One")
loadFiles <- list.files(path = "/nfs/turbo/umms-scjp/kgmoo/data/rna_post_decontx_counts", pattern = "*.Rda", full.names = TRUE, recursive = TRUE)
groupFiles <- data.frame(t(data.frame(strsplit(loadFiles, split = "/"))))$X8
groupNames <- data.frame(t(data.frame(strsplit(groupFiles, split = ".d"))))$X1
geneList <- read.table(file = "/nfs/turbo/umms-scjp/kgmoo/data/rna_post_decontx_counts/692-NM-3.decont/genes.tsv", sep = "\t")

# Arushi indicated this file will give us the nuclei whitelist
print("Two")
metadataFactors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/cluster_info_qc.tsv", sep = "\t", header = TRUE)

# Include only rna nuclei (not atac and multiomic nuclei)
metadataFactors <- metadataFactors[(metadataFactors$modality == "rna") & (metadataFactors$demux_cat != "multiome"),]

# First fence post loaded after removing libraries that should not be included
print("Three")
loadFiles <- loadFiles[-1] # Hard code drop the multiome file from the list, perhaps not the best solution
load(loadFiles[1])
sceDHEAVY <- m1
dim(sceDHEAVY)
loadFiles <- loadFiles[-1]
sceDONE <- sceDHEAVY

# Finish appending rna counts files together
print("Four")
for (tempFile in loadFiles) {
    load(file = tempFile)
    tempSCE <- m1
    sceDHEAVY <- cbind(sceDHEAVY, tempSCE, deparse.level=1)    
}

# Generate the cell level unique id column that Arushi's metadata files tend to use to index the results
print("Five")
sceDHEAVY$id <- paste0("rna.",gsub("-", ".", sceDHEAVY$Sample),".",sceDHEAVY$Barcode)

sceDHEAVY_Filter <- sceDHEAVY[,sceDHEAVY$id %in% metadataFactors$index]

# Load the more complete metadata file so we can get the model into place, select the rows and columns useful to us
print("Six")
sampMet <- read.table(file = "/nfs/turbo/umms-scjp/kgmoo/data/fusion_sample_info.tsv", header = TRUE, sep = "\t")

# Build complete metadata frame
print("Seven")
colOut <- colData(sceDHEAVY)
colOut <- data.frame(merge(colOut, metadataFactors, by.x = "id", by.y = "index")) %>% 
dplyr::select(c(id, Barcode, UMAP_1, UMAP_2, modality, library, SNG.1ST, hqaa_umi, fraction_mitochondrial, sample, coarse_cluster_name))
colOut <- data.frame(merge(colOut, sampMet, by = "SNG.1ST"))
colOut <- colOut %>% dplyr::arrange(id)

# Attach desired metadata to the sce object
sceDHEAVY_Filter$UMAP_1 <- as.numeric(colOut$UMAP_1)
sceDHEAVY_Filter$UMAP_2 <- as.numeric(colOut$UMAP_2)
sceDHEAVY_Filter$modality <- as.factor(colOut$modality)
sceDHEAVY_Filter$hqaa_umi <- colOut$hqaa_umi
sceDHEAVY_Filter$fraction_mitochondrial <- colOut$fraction_mitochondrial
sceDHEAVY_Filter$sample <- colOut$sample
sceDHEAVY_Filter$coarse_cluster_name <- colOut$coarse_cluster_name
sceDHEAVY_Filter$age <- as.numeric(colOut$age)
sceDHEAVY_Filter$sex <- as.factor(colOut$sex)
sceDHEAVY_Filter$batch <- as.factor(colOut$batch)
sceDHEAVY_Filter$bmi <- as.numeric(colOut$bmi)

# For each cell type, aggregate counts per gene by cell to get a total counts per sample sheet like the one used for peak tracking
clustersOI <- c('Type_1','Type_2x','Type_2a','Endothelial','Mesenchymal_Stem_Cell','Macrophage','Satellite_Cell','Neuromuscular_junction','Smooth_Muscle','Neuronal','Muscle_Fiber_Mixed','T_cell','Adipocyte')

sampleList <- unique(sceDHEAVY_Filter$sample)
for (cDex in clustersOI) {
    print(cDex)
    tempRes <- data.frame(Ens = rowData(sceDHEAVY_Filter)$ID, Symbol = rowData(sceDHEAVY_Filter)$Symbol)
    
    # Processing the cell type specific SCE for seperate export as SCA for MAST
    # We filter for lowly expressed genes and samples with few to no cells
    print("Eight")
    sceDTYPE <- sceDHEAVY_Filter[,sceDHEAVY_Filter$coarse_cluster_name == cDex]
    countsTemp <- counts(sceDTYPE)
    sampleTemp <- colData(sceDTYPE)$sample
    sampleTempList <- unique(sampleTemp)

    # Calculate CDR as recommended by MAST documentation 
    sceDTYPE$CDR <- as.numeric(colSums(counts(sceDTYPE)>0)/nrow(counts(sceDTYPE)))
    print(min(sceDTYPE$CDR))
    print(max(sceDTYPE$CDR))

    # This filter is not used as currently written, but provides an alternative filter done on the raw counts
    filterFeat <- rowSums(countsTemp >= 1) >= ncol(countsTemp)*0.1 #length(unique(colData(sceDTYPE)$sample))*0.1
    sceDTYPEfor <- sceDTYPE#[filterFeat,] #TPM moved to later upon request
    print("Nine")
      
    # Use decontXcounts as the basis and then calculate log2(TPM+1) counts as required by MAST. Finally apply filter from earlier
    counts(sceDTYPEfor) <- assay(sceDTYPEfor, "decontXcounts")
    counts(sceDTYPEfor) <- calculateTPM(sceDTYPEfor)
    logcounts(sceDTYPEfor) <- log2(counts(sceDTYPEfor)+1)
    scaDTYPE <- SceToSingleCellAssay(sceDTYPEfor, check_sanity = FALSE)
    scaDTYPEF <- subset(scaDTYPE, sex == "F")
    scaDTYPEM <- subset(scaDTYPE, sex == "M")
    gp <- freq(scaDTYPE)
    gpF <- freq(scaDTYPEF)
    gpM <- freq(scaDTYPEM)

    backup <- scaDTYPE
    backupF <- scaDTYPEF
    backupM <- scaDTYPEM
    
    tryCatch({
            genes_to_keep <- rownames(backup[gp >= 0.1,])   
        }, error = function(e){
            print(e)
            genes_to_keep <- c()
        })
    tryCatch({
            genes_to_keepF <- rownames(backupF[gpF >= 0.1,])   
        }, error = function(e){
            print(e)
            genes_to_keepF <- c()
        })
    tryCatch({
            genes_to_keepM <- rownames(backupM[gpM >= 0.1,])   
        }, error = function(e){
            print(e)
            genes_to_keepM <- c()
        })
    
    
    scaDTYPE <- scaDTYPE[rownames(scaDTYPE) %in% genes_to_keep,]
    scaDTYPEF <- scaDTYPEF[rownames(scaDTYPEF) %in% genes_to_keepF,]
    scaDTYPEM <- scaDTYPEM[rownames(scaDTYPEM) %in% genes_to_keepM,]
    #scaDTYPEfor <- scaDTYPE[filterFeat,]

    saveRDS(scaDTYPE, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nTotal_CDR_sca.rds"))
    saveRDS(scaDTYPEF, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nFemale_CDR_sca.rds"))
    saveRDS(scaDTYPEM, file = glue("/nfs/turbo/umms-scjp/kgmoo/data/rna_sample_counts/{cDex}_rna_postDecontx_postArushiQC_PostTPM_1in10nMale_CDR_sca.rds"))
    
}

print("Done")
