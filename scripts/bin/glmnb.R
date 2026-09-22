#! /usr/bin/env Rscript
suppressPackageStartupMessages({library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(glue)
library(cowplot)
library(pbmcapply)
library(DHARMa)
library(magrittr)
library(glmmTMB)
library(pscl)
library(optparse)
library(lmtest)
})

option_list <- list(
    make_option(c("--cluster"), type = "character", help = "[Required] cluster to select"),
    make_option(c("--contrast"), type = "character", help = "[Required] colon-separated contrast  eg No-HFPEF:HFPEF"),
    make_option(c("--sce"), type = "character", help = "[Required] ATAC single cell experiment object"),
    make_option(c("--n_threads"), type = "numeric", default = 1, help = "Number of cores to use"),
    make_option(c("--n_chunks"), type = "numeric", default = 1, help = "Number of chunk to split peaks into while running the model"),
    make_option(c("--chunk_number"), type = "numeric", default = 1, help = "nth chunk to run on"),
    make_option(c("--prefix"), type = "character", help = "[Required] output prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# Establish a function to save images for the script
savep = function(p, filename, height=6, width=7){
  png(filename = filename, width = width, height = height, units = "in", res = 150 )
  print(p)
  dev.off()
  print(p)
}

# Load in all peaks and subset to the cluster currently being analyzed 
#cluster = opts$cluster
#merged_peaks= read.table(opts$peaks, sep='\t', as.is=T, header=T)
#selected.peaks = merged_peaks[merged_peaks[[cluster]]=="true",]$Peaks
### Counts or just coords??? ###

small <- readRDS(opts$sce)

#small$CDR <- as.numeric(colSums(counts(small)>0)/nrow(counts(small)))

smallF <- small[,small$sex == "F"]
smallM <- small[,small$sex == "M"]
contrast <- unlist(strsplit(opts$contrast, ":"))
### test with default sce method #AV#

# Apply filters to genes for testing
keepRows <- rowSums(counts(small) >= 1) >= ncol(counts(small))*0.1
small <- small[keepRows,]
keepRowsF <- rowSums(counts(smallF) >= 1) >= ncol(counts(smallF))*0.1
smallF <- smallF[keepRowsF,]
keepRowsM <- rowSums(counts(smallM) >= 1) >= ncol(counts(smallM))*0.1
smallM <- smallM[keepRowsM,]

# Extract the count data from the sce and drop all of the factors that are now no longer represented
expr <- counts(small)
exprF <- counts(smallF)
exprM <- counts(smallM)

meta <- droplevels(colData(small))
meta$batch = as.factor(meta$batch)
meta$ogtt = as.factor(meta$ogtt)
meta$bmi = as.numeric(meta$bmi)
meta$age = as.numeric(meta$age)
meta$tss_enrichment = as.numeric(meta$tss_enrichment)
#meta$CDR = as.numeric(meta$CDR)

metaF <- meta[meta$sex == "F",]
metaM <- meta[meta$sex == "M",]

meta$sex = as.factor(meta$sex)

# Set test and null formula where the only difference is the covariate being tested
# This allows your to compare the two via likelhood raitio test to get a p-value
if (opts$contrast %in% c("No-HFPEF:HFPEF", "No-HFPEF:HFPEF-followup")){
  fmla <- as.formula("GENE ~ group + age + batch + sex_label + offset(log(total_counts)) + (1|samp)")
  null_fmla <- as.formula("GENE ~ age + batch + sex_label + offset(log(total_counts)) + (1|samp)")
} else if (opts$contrast == "HFPEF:HFPEF-followup"){
  fmla <- as.formula("GENE ~ group + batch + individual + offset(log(total_counts)) + (1|samp)")
  null_fmla <- as.formula("GENE ~ batch + individual + offset(log(total_counts)) + (1|samp)")
} else if (opts$contrast == "age"){
  print("Age Contrast IDd")
  fmla <- as.formula("GENE ~ age + bmi + batch + sex + ogtt + tss_enrichment + offset(log(total_counts)) + (1|individual)")
  null_fmla <- as.formula("GENE ~ bmi + batch + sex + ogtt + tss_enrichment + offset(log(total_counts))+ (1|individual)")
  fmla_SS <- as.formula("GENE ~ age + bmi + batch + ogtt + tss_enrichment + offset(log(total_counts))+ (1|individual)")
  null_fmla_SS <- as.formula("GENE ~ bmi + batch + ogtt + tss_enrichment + offset(log(total_counts))+ (1|individual)")
} else {
  print("Contrast Not Found")
}

# Build a function for running the model
glmnb = function(index, expr, meta, fmla, null_fmla){
  # Takes the expression from one gene at a time and blocks it together with metadata
  # I will need to add ogtt and bmi here after ensuring they are part of meta above
  print("index")
  print(index)
  test_dat = data.frame(
    GENE = expr[index,],
    samp = meta$sample,
    batch = as.factor(meta$batch),
    sex = as.factor(meta$sex),
    ogtt = as.factor(meta$ogtt),
    bmi = as.numeric(meta$bmi),
    age = as.numeric(meta$age),
    tss_enrichment = as.numeric(meta$tss_enrichment)
    #CDR = as.numeric(meta$CDR)
  )
  # Generate new columns for the model, total_counts of this one gene (+1?)
  # and an individual id that includes the batch number... it seems. 
  test_dat %<>% mutate(total_counts = as.numeric(colSums(expr)+1))
  test_dat %<>% mutate(individual = mapply(gsub, glue("--{batch}"), "", samp))
  
  ## initialize with NA in case of error late #AV#
  age = estimate = se = z = p_mod1 = p_val = NA
  print("Begin Model ---------------------------------------------------")
  #print(fmla) # -----------------------------------------------------------
  #print(null_fmla) # ------------------------------------------------------
  #print(head(test_dat)) # -------------------------------------------------

  tryCatch({
    print("Try Catch Entered ----------------------------------------------")
    # Run the test and null models using the glmmTMB method from package
    mod1 = glmmTMB(fmla, test_dat, family = nbinom2, REML = FALSE)
    mod_null = glmmTMB(null_fmla, test_dat, family = nbinom2, REML = FALSE)
    #print(head(summary(mod1)$coefficients$cond)) # ------------------------
    # Get the gene names, coeffecients, standard errors, z values, and pvalues of the age component of the model
    hfpf_index = grep("age", rownames(summary(mod1)$coefficients$cond))
    hfpf = rownames(summary(mod1)$coefficients$cond)[hfpf_index]
    estimate = summary(mod1)$coefficients$cond[hfpf,'Estimate']
    se = summary(mod1)$coefficients$cond[hfpf,'Std. Error']
    z = summary(mod1)$coefficients$cond[hfpf,'z value']
    p_mod1 = summary(mod1)$coefficients$cond[hfpf,'Pr(>|z|)']

    # Run the likelihood ratio test comparing the test and null models, get that p-value
    lr_result = lrtest(mod1, mod_null)
    p_val <- lr_result$`Pr(>Chisq)`[2] 
    print(p_val)

  }, error = function(e) {print(e)})
  result = list(index, hfpf, estimate, se, z, p_mod1, p_val) 
  print(result)
  return(result)
  
}


# Get number of chunks and assign an array of them to work on each peak
n_chunks = opts$n_chunks
chunk_number = opts$chunk_number
npeaks = length(rownames(expr))
npeaksF <- length(rownames(exprF))
npeaksM <- length(rownames(exprM))

len_chunk = round(npeaks/n_chunks)
len_chunkF <- round(npeaksF/n_chunks)
len_chunkM <- round(npeaksM/n_chunks)
start = (chunk_number - 1) * (len_chunk) + 1
startF = (chunk_number - 1) * (len_chunkF) + 1
startM = (chunk_number - 1) * (len_chunkM) + 1
end = start + len_chunk
endF = startF + len_chunkF
endM = startM + len_chunkM
if (chunk_number == n_chunks){
  end = npeaks
}
if (chunk_number == n_chunks){
  endF = npeaksF
}
if (chunk_number == n_chunks){
  endM = npeaksM
}

# Set seed and run all peaks that pass the filter
set.seed(1055)
#print("Combined ********************************")
peaks_to_run = rownames(expr)[start:end]
print(glue("selecting {n_chunks} out of {npeaks}"))
 #When testing skip the pbmclapply as that is for the multithreading
 #p_GLMMNB is just for your contrst age
 #p_LRT is the p-value comparing the model vs the null model, you want this one for downstream use
 #You may need to filter out NA results and the corresponding peaks
out = pbmclapply(peaks_to_run, mc.cores = opts$n_threads, function(i){glmnb(i, expr, meta, fmla, null_fmla)})
res = do.call(rbind, lapply(out, function(i){as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", "p_LRT"))}))
#res = do.call(rbind, lapply(out, function(i){as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", 
#"estimate_sex", "se_sex", "z_sex", "p_GLMMNB_sex", "p_LRT"))}))

#print("Female ***********************************")
#peaks_to_runF = rownames(exprF)[startF:endF]
#print(glue("selecting {n_chunks} out of {npeaksF}"))
#outF = pbmclapply(peaks_to_runF, mc.cores = opts$n_threads, function(i){glmnb(i, exprF, metaF, fmla_SS, null_fmla_SS)})
#resF = do.call(rbind, lapply(outF, function(i){as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", "p_LRT"))}))

#print("Male *************************************")
#peaks_to_runM = rownames(exprM)[startM:endM]
#print(glue("selecting {n_chunks} out of {npeaksM}"))
#outM = pbmclapply(peaks_to_runM, mc.cores = opts$n_threads, function(i){glmnb(i, exprM, metaM, fmla_SS, null_fmla_SS)})
#resM = do.call(rbind, lapply(outM, function(i){as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", "p_LRT"))}))

#print(head(resF))
write.table(res, sep='\t', quote = FALSE, file = glue("pTSS_{opts$prefix}.tsv"))
#write.table(resF, sep='\t', quote = FALSE, file = glue("F_pTSS_{opts$prefix}.tsv"))
#write.table(resM, sep='\t', quote = FALSE, file = glue("M_pTSS_{opts$prefix}.tsv"))
