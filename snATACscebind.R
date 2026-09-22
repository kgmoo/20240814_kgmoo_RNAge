##################
### 2025/05/19 ###
### ATAC Cast  ###
##################

# This script loads in long form libraries of snATAC-seq results and reads them into an sce object.
# It is incredibly memory intensive as it must shift the data between two different memory efficient
# sparse datatype storage. In this case no subsetting was done as Steve asked for all nuclei type
# to eventually be run. 

suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("data.table"))

cDex1 <- "Neuromuscular_junction"
cDex2 <- "Satellite_Cell"
cDex3 <- "Muscle_Fiber_Mixed"
cDex4 <- "T_cell"

# Remove the multiome library
libDir <- "/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/"
libList <- list.files("/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/")[-1]
libInit <- libList[1]
libList <- libList[-1]

allPeaks <- readRDS("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/allPeaks.rds")

metadata.factors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/cluster_info_qc.tsv", header = TRUE) %>% filter(modality == "atac")
print("metadata")
print(dim(metadata.factors))
print(head(metadata.factors))
print(unique(metadata.factors$coarse_cluster_name))
sampledata.factors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/batch-sample-info-full.tsv", sep = "\t", header = TRUE) %>% select(SNG.1ST, sex, age, bmi, ogtt)
print(head(sampledata.factors))
metadata.factorSamp <- merge(metadata.factors, sampledata.factors, by = "SNG.1ST")

seedFrame <- data.frame(read.table(file = glue(libDir, libInit), sep="\t"))
colnames(seedFrame) <- c("peakV", "nucV", "countV")
initNucs <- sort(unique(seedFrame$nucV))
 
metadata.factorSamp_seed <- merge(data.frame(index = initNucs), metadata.factorSamp)

seedSparse <- sparseMatrix(i=match(seedFrame$peakV, allPeaks), j=match(seedFrame$nucV, initNucs), x = seedFrame$countV)
print("seedSparse")
print(dim(seedSparse))
rownames(seedSparse) <- allPeaks
colnames(seedSparse) <- initNucs

sce <- SingleCellExperiment(list(counts=seedSparse), colData=metadata.factorSamp_seed)
sce1 <- sce[,sce$coarse_cluster_name == cDex1]
sce2 <- sce[,sce$coarse_cluster_name == cDex2]
sce3 <- sce[,sce$coarse_cluster_name == cDex3]
sce4 <- sce[,sce$coarse_cluster_name == cDex4]
print("seedSce")
print(dim(sce1))
for (file in libList) {
    
    tempFrame <- data.frame(read.table(file = glue(libDir, file), sep="\t"))
    colnames(tempFrame) <- c("peakV", "nucV", "countV")
    tempNucs <- sort(unique(tempFrame$nucV))
    
    print(file)
    print(dim(tempFrame))
    
    metadata.factorSamp_temp <- merge(data.frame(index = tempNucs), metadata.factorSamp)
    
    tempSparse <- sparseMatrix(i=match(tempFrame$peakV, allPeaks), j=match(tempFrame$nucV, tempNucs), x = tempFrame$countV)
    print("tempSparse")
    print(dim(tempSparse))
    rownames(tempSparse) <- allPeaks
    colnames(tempSparse) <- tempNucs
    
    tempSce <- SingleCellExperiment(list(counts=tempSparse), colData=metadata.factorSamp_temp)
    tempSce1 <- tempSce[,tempSce$coarse_cluster_name == cDex1]
    tempSce2 <- tempSce[,tempSce$coarse_cluster_name == cDex2]
    tempSce3 <- tempSce[,tempSce$coarse_cluster_name == cDex3]
    tempSce4 <- tempSce[,tempSce$coarse_cluster_name == cDex4]
    print("tempSce")
    print(dim(tempSce1))
    
    sce1 <- cbind(sce1, tempSce1)
    sce2 <- cbind(sce2, tempSce2)
    sce3 <- cbind(sce3, tempSce3)
    sce4 <- cbind(sce3, tempSce4)
    print("new sce")
    print(dim(sce1))
    
}

print("Done")
print(dim(sce1))

saveRDS(sce1, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/{cDex1}_FUSIONsnATAC_sce.rds"))
saveRDS(sce2, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/{cDex2}_FUSIONsnATAC_sce.rds"))
saveRDS(sce3, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/{cDex3}_FUSIONsnATAC_sce.rds"))
saveRDS(sce4, glue("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/{cDex4}_FUSIONsnATAC_sce.rds"))
