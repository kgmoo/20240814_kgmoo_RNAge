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
    make_option(c("--peaks"), type = "character", help = "[Required] dataframe with peak calls in clusters"),
    make_option(c("--contrast"), type = "character", help = "[Required] colon-separated contrast  eg No-HFPEF:HFPEF"),
    make_option(c("--sce"), type = "character", help = "[Required] ATAC single cell experiment object"),
    make_option(c("--n_threads"), type = "numeric", default = 1, help = "Number of cores to use"),
    make_option(c("--n_chunks"), type = "numeric", default = 1, help = "Number of chunk to split peaks into while running the model"),
    make_option(c("--chunk_number"), type = "numeric", default = 1, help = "nth chunk to run on"),
    make_option(c("--prefix"), type = "character", help = "[Required] output prefix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

savep = function(p, filename, height=6, width=7){
  png(filename = filename, width = width, height = height, units = "in", res = 150 )
  print(p)
  dev.off()
  print(p)
}

cluster = opts$cluster
merged_peaks= read.table(opts$peaks, sep='\t', as.is=T, header=T)
selected.peaks = merged_peaks[merged_peaks[[cluster]]=="true",]$Peaks

load(opts$sce)
contrast = unlist(strsplit(opts$contrast, ":"))
### test with default sce method
small = pads[selected.peaks, colData(pads)$group %in% contrast & colData(pads)$coarse.clust %in% c(cluster)]
if (opts$contrast == "No-HFPEF:HFPEF"){
  small = small[,!grepl("6721-NM-56", colData(small)$samp)]
  print("Removed 6721-NM-56 as it's the same individual as 6721-NM-17")
}

expr = assay(small)
meta =  droplevels(colData(small))
meta$batch = as.factor(meta$batch)
meta$sex_label = as.factor(meta$sex_label)
print(head(meta))

if (opts$contrast %in% c("No-HFPEF:HFPEF", "No-HFPEF:HFPEF-followup")){
  fmla <- as.formula("GENE ~ group + age + batch + sex_label + offset(log(total_counts)) + (1|samp)")
  null_fmla <- as.formula("GENE ~ age + batch + sex_label + offset(log(total_counts)) + (1|samp)")
  fmla_random <- as.formula("GENE ~ group + age + (1|batch) + sex_label + offset(log(total_counts)) + (1|samp)")
} else if (opts$contrast == "HFPEF:HFPEF-followup"){
  fmla <- as.formula("GENE ~ group + batch + individual + offset(log(total_counts)) + (1|samp)")
  null_fmla <- as.formula("GENE ~ batch + individual + offset(log(total_counts)) + (1|samp)")
  fmla_random <- as.formula("GENE ~ group + batch + (1|individual) + offset(log(total_counts)) + (1|samp)")
} 

glmnb = function(index, fmla, null_fmla, fmla_random){
  test_dat = data.frame(
    GENE = expr[index,],
    group = meta$group,
    samp = meta$samp,
    batch = as.factor(meta$batch),
    sex_label = as.factor(meta$sex_label),
    age = meta$age
  )
  test_dat %<>% mutate(total_counts = as.numeric(colSums(expr)+1))
  # print(head(test_dat))
  test_dat %<>% mutate(individual = mapply(gsub, glue("-{batch}"), "", samp))
  # print(head(test_dat))
  tryCatch({
    mod1 = glmmTMB(fmla, test_dat, family = nbinom2, REML = FALSE)
    mod_null = glmmTMB(null_fmla, test_dat, family = nbinom2, REML = FALSE)
    mod_random = glmmTMB(fmla_random, test_dat, family = nbinom2, REML = FALSE)

    hfpf_index = grep("HFPEF", rownames(summary(mod1)$coefficients$cond))
    hfpf = rownames(summary(mod1)$coefficients$cond)[hfpf_index]
    estimate = summary(mod1)$coefficients$cond[hfpf,'Estimate']
    se = summary(mod1)$coefficients$cond[hfpf,'Std. Error']
    z = summary(mod1)$coefficients$cond[hfpf,'z value']
    p_mod1 = summary(mod1)$coefficients$cond[hfpf,'Pr(>|z|)']

    # estimate_sex = summary(mod1)$coefficients$cond['sex_labelMale','Estimate']
    # se_sex = summary(mod1)$coefficients$cond['sex_labelMale','Std. Error']
    # z_sex = summary(mod1)$coefficients$cond['sex_labelMale','z value']
    # p_mod1_sex = summary(mod1)$coefficients$cond['sex_labelMale','Pr(>|z|)']
    print(summary(mod1))

    lr_result = lrtest(mod1, mod_null)
    p_val <- lr_result$`Pr(>Chisq)`[2] 

    lr_result = lrtest(mod1, mod_random)
    p_val_random <- lr_result$`Pr(>Chisq)`[2] 
    
    sim_res <- simulateResiduals(mod1)
    p_zero = testZeroInflation(sim_res)$p.value

    result = list(index, hfpf, estimate, se, z, p_mod1, p_val, p_val_random, p_zero) 
    # result = list(index, hfpf, estimate, se, z, p_mod1, estimate_sex, se_sex, z_sex, p_mod1_sex, p_val)

  }, error = function(e) {print(e); result = list(index, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)})
  return(result)
  
}

n_chunks = opts$n_chunks
# chunk_number = opts$chunk_number
npeaks = length(rownames(expr))
# len_chunk = round(npeaks/n_chunks)
# start = (chunk_number - 1) * (len_chunk) + 1
# end = start + len_chunk
# if (chunk_number == n_chunks){
#   end = npeaks
# }

set.seed(1234)
peaks_to_run = sample(rownames(expr), n_chunks)
# peaks_to_run = rownames(expr)[start:end]
print(glue("selecting {n_chunks} out of {npeaks}"))
out = pbmclapply(peaks_to_run, mc.cores = opts$n_threads, function(i){glmnb(i, fmla, null_fmla, fmla_random)})
res = do.call(rbind, lapply(out, function(i){as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", "p_LRT", "p_random", "p_zero"))}))
# res = do.call(rbind, lapply(out, function(i){as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", 
# "estimate_sex", "se_sex", "z_sex", "p_GLMMNB_sex", "p_LRT"))}))

print(res)
write.table(res, sep='\t', quote = FALSE, file = glue("{opts$prefix}.tsv"))

## fold change
sc0 = CreateSeuratObject(expr, meta.data = meta) %>%  NormalizeData()
mat = GetAssayData(sc0, slot = 'data')
if (opts$contrast == "No-HFPEF:HFPEF"){
  cells1 = rownames(meta[meta$group=="No-HFPEF",])
  cells2 = rownames(meta[meta$group=="HFPEF",])
} else if (opts$contrast == "HFPEF:HFPEF-followup"){
  cells1 = rownames(meta[meta$group=="HFPEF",])
  cells2 = rownames(meta[meta$group=="HFPEF-followup",]) 
} else if (opts$contrast == "No-HFPEF:HFPEF-followup"){
  cells1 = rownames(meta[meta$group=="No-HFPEF",])
  cells2 = rownames(meta[meta$group=="HFPEF-followup",]) 
}
# cells1 = rownames(meta[meta$group=="No-HFPEF",])
# cells2 = rownames(meta[meta$group=="HFPEF",])
data1 = log(rowMeans(expm1(mat[, cells1, drop = F]) + 1))
data2 = log(rowMeans(expm1(mat[, cells2, drop = F]) + 1))
lfc = as.numeric(data2 - data1)
lfc = data.frame(lfc)
lfc$index = rownames(mat)

res = merge(res, lfc, by="index")

write.table(res, sep='\t', quote = FALSE, file = glue("{opts$prefix}.tsv"))