##################
### Keagan Moo ###
### 2025/05/19 ###
### ATAC Cast  ###
##################

## DEPRECATED - Use snATACscebind.R

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

# Remove the multiome library
libDir <- "/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/"
libList <- list.files("/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/")[-1]
libInit <- libList[1]
libList <- libList[-1]

allPeaks <- readRDS("/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/allPeaks.rds")

seedFrame <- data.frame(read.table(file = glue(libDir, libInit), sep="\t"))
colnames(seedFrame) <- c("peakV", "nucV", "countV")
initNucs <- sort(unique(seedFrame$nucV))

seedSparse <- sparseMatrix(i=match(seedFrame$peakV, allPeaks), j=match(seedFrame$nucV, initNucs), x = seedFrame$countV)
print("seedSparse")
print(dim(seedSparse))
rownames(seedSparse) <- allPeaks
colnames(seedSparse) <- initNucs

for (file in libList) {
    
    tempFrame <- data.frame(read.table(file = glue(libDir, file), sep="\t"))
    colnames(tempFrame) <- c("peakV", "nucV", "countV")
    tempNucs <- sort(unique(tempFrame$nucV))
    
    print(file)
    print(dim(tempFrame))
    
    tempSparse <- sparseMatrix(i=match(tempFrame$peakV, allPeaks), j=match(tempFrame$nucV, tempNucs), x = tempFrame$countV)
    print("tempSparse")
    print(dim(tempSparse))
    rownames(tempSparse) <- allPeaks
    colnames(tempSparse) <- tempNucs
    
    seedSparse <- cbind(seedSparse, tempSparse)
    print("new seedSparse")
    print(dim(seedSparse))
    
}

print("Done")
print(dim(seedSparse))
print(head(seedSparse))

metadata.factors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/cluster_info_qc.tsv", header = TRUE) %>% filter(modality == "atac")
print("metadata")
print(dim(metadata.factors))
print(head(metadata.factors))
sampledata.factors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/batch-sample-info-full.tsv", sep = "\t", header = TRUE) %>% select(SNG.1ST, sex, age, bmi, ogtt)
print(head(sampledata.factors))
metadata.factorSamp <- merge(metadata.factors, sampledata.factors, by = "SNG.1ST")

allNucs <- sort(unique(colnames(seedSparse)))
 
metadata.factorSamp <- merge(data.frame(index = allNucs), metadata.factorSamp)

sce <- SingleCellExperiment(list(counts=seedSparse), colData=metadata.factorSamp)

saveRDS(sce, "/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/FUSIONsn_sce.rds")
