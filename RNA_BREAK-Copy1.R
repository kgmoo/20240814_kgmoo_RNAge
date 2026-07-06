##################
### Keagan Moo ###
### 2024/08/15 ###
### RNA Break  ###
##################

# This script loads in the nuclei and sample level data from the single nucleus skeletal muscle data from FUSION,
# and using the MILO environment breaks them up into cell type cell x gene matricies for later work
library("Seurat")
library("glue")
library("SingleCellExperiment")
library("scater")
library("scran")
library("tidyverse")
library("data.table")
library("MAST")

# Get list of rna count files as well as genes that occur in them
print("One")
loadFiles <- list.files(path = "/nfs/turbo/umms-scjp/kgmoo/data/rna_post_decontx_counts", pattern = "*.Rda", full.names = TRUE, recursive = TRUE)
groupFiles <- data.frame(t(data.frame(strsplit(loadFiles, split = "/"))))$X8
groupNames <- data.frame(t(data.frame(strsplit(groupFiles, split = ".d"))))$X1
geneList <- read.table(file = "/nfs/turbo/umms-scjp/kgmoo/data/rna_post_decontx_counts/692-NM-3.decont/genes.tsv", sep = "\t")

# Arushi indicated this file will give us the nuclei that should be included in the analysis after qc
print("Two")
metadataFactors <- read.table("/nfs/turbo/umms-scjp/kgmoo/data/cluster_info_qc.tsv", sep = "\t", header = TRUE)
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

# Append desired metadata to the sce object
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

print(colnames(colData(sceDHEAVY_Filter)))
forCollapse <- data.frame(colData(sceDHEAVY_Filter))

print("UMI per nuclei")
mean(forCollapse$hqaa_umi)
median(forCollapse$hqaa_umi)

print("Genes per nuclei")
mean(colSums(counts(sceDHEAVY_Filter) > 0))
median(colSums(counts(sceDHEAVY_Filter) > 0))

readsDF <- data.frame(names(colSums(counts(sceDHEAVY_Filter))), colSums(counts(sceDHEAVY_Filter)))
colnames(readsDF) <- c("Barcode", "nReads")

forCollapse <- merge(forCollapse, readsDF, by = c("Barcode"))

sampledf <- forCollapse %>%
  group_by(Sample) %>%
  summarise(reads = sum(nReads), .groups = "drop")
                     
print("Reads per sample")
mean(sampledf$reads)
median(sampledf$reads)

forCollapse <- forCollapse[forCollapse$coarse_annot == "Type_2a",]
sampledf <- forCollapse %>%
group_by(Sample) %>%
summarise(reads = sum(nReads), .groups = "drop")
                     
print("Reads per sample (Cell Type)")
mean(sampledf$reads)
median(sampledf$reads)