##################
### Keagan Moo ###
### 2025/05/21 ###
### ATAC Cast  ###
##################

# This script simply gets all unique peaks from the different ATAC libraries and compiles them into a single sorted list
# so that sparse matricies can be generated on a per library basis with matching rows

suppressPackageStartupMessages({
library(Matrix)
library(dplyr)
library(glue)
})

# Load all libraries, remove the multiome library, as well as the first library which will be used as a seed vector for the loop

libDir <- "/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/"
libList <- list.files("/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/")[-1]
libInit <- libList[1]
libList <- libList[-1]

# Load the initial libarary's unqiue peak list and order it
seedLib <- sort(unique(read.table(
    '/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/quantify_atac/results/fragment-counts/fusion-150/499-NM-1.counts.txt', sep = "\t")$V1))

print("499-NM-1")
print(length(seedLib))
# For each remaining library
for (file in libList) {
    print(file)
    
    # Load just the unique sorted peak list from the library
    tempLib <- sort(unique(read.table(
    glue(libDir, file), sep = "\t")$V1))
    print(length(tempLib))
    
    # Cat the complete vector of peaks with the new unique list and reassert unique members and sorting in that order
    seedLib <- sort(unique(c(seedLib, tempLib)))
    print("Merged...")
    print(length(seedLib))
}

print(head(seedLib))
# Save the vector for loading into the full sce assembly script
saveRDS(seedLib, file = "/nfs/turbo/umms-scjp/kgmoo/analysis/20250428_kgmoo_ATACge/data/allPeaks.rds")