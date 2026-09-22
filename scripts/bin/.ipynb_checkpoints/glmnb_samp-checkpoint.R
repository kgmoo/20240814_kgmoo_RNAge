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
    make_option(c("--type"), type = "character", default="coarse.clust", help = "[Required] type of peaks"),
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
type = opts$type

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
print(head(colData(pads)[[type]]))
#small = pads[selected.peaks, colData(pads)$group %in% contrast & colData(pads)[[type]] == cluster]

expr = assay(small)
meta =  droplevels(colData(small))
meta$batch = as.factor(meta$batch)
meta$sex_label = as.factor(meta$sex_label)
print(head(meta))


fmla <- as.formula("GENE ~ group + batch + individual + offset(log(total_counts)) + (1|samp)")
null_fmla <- as.formula("GENE ~ batch + individual + offset(log(total_counts)) + (1|samp)")

fmla_random <- as.formula("GENE ~ group + batch + (1|individual) + offset(log(total_counts)) + (1|samp)")
null_fmla_random <- as.formula("GENE ~ batch + (1|individual) + offset(log(total_counts)) + (1|samp)")

fitglm = function(fmla, null_fmla, test_dat){
  ## initialize with NA in case of error late
  hfpf = estimate = se = z = p_mod = p_lrt = NA

  mod = glmmTMB(fmla, test_dat, family = nbinom2, REML = FALSE)
  mod_null = glmmTMB(null_fmla, test_dat, family = nbinom2, REML = FALSE)

  hfpf_index = grep("HFPEF", rownames(summary(mod)$coefficients$cond))
  hfpf = rownames(summary(mod)$coefficients$cond)[hfpf_index]
  estimate = summary(mod)$coefficients$cond[hfpf,'Estimate']
  se = summary(mod)$coefficients$cond[hfpf,'Std. Error']
  z = summary(mod)$coefficients$cond[hfpf,'z value']
  p_mod = summary(mod)$coefficients$cond[hfpf,'Pr(>|z|)']

  print(summary(mod))

  lr_result = lrtest(mod, mod_null)
  p_lrt <- lr_result$`Pr(>Chisq)`[2] 
  return(list(hfpf = hfpf, estimate = estimate, se = se, z = z, p_mod = p_mod, p_lrt = p_lrt))
}

glmnb = function(index, fmla, null_fmla, fmla_random, null_fmla_random){
  test_dat = data.frame(
    GENE = expr[index,],
    group = meta$group,
    samp = meta$samp,
    batch = as.factor(meta$batch),
    sex_label = as.factor(meta$sex_label),
    study = as.factor(meta$study),
    age = meta$age
  )
  test_dat %<>% mutate(total_counts = as.numeric(colSums(expr)+1))
  test_dat %<>% mutate(individual = mapply(gsub, glue("-{batch}"), "", samp))
  # head(test_dat)
  tryCatch({
    res = fitglm(fmla_random, null_fmla_random, test_dat)
    ncontrol = length(unique(test_dat[test_dat$group=="HFPEF",]$individual))
    ncase = length(unique(test_dat[test_dat$group=="HFPEF-followup",]$individual))
    result_random = list(index, res$hfpf, res$estimate, res$se, res$z, res$p_mod, res$p_lrt, "random_individual", ncase, ncontrol)
    
    ## subset to paired only samples and run fixed model
    test_dat = test_dat[test_dat$study == "SECRETII",]
    ncontrol = length(unique(test_dat[test_dat$group=="HFPEF",]$individual))
    ncase = length(unique(test_dat[test_dat$group=="HFPEF-followup",]$individual))
    res1 = fitglm(fmla, null_fmla, test_dat)
    result_fixed = list(index, res1$hfpf, res1$estimate, res1$se, res1$z, res1$p_mod, res1$p_lrt, "fixed_individual_paired", ncontrol, ncase)
    }, error = function(e) {print(e)})
  # print(result_random)
  print(result_fixed)

  output <- do.call(rbind, lapply(list(result_random, result_fixed), function(i){
    as.data.frame(i, col.names = c("index", "group_i", "estimate", "se", "z", "p_GLMMNB", "p_LRT", "model", "nsamp_hfpef", "nsamp_hfpef-followup"))
    }))
  print(output)
  return(output)
  
}

n_chunks = opts$n_chunks
chunk_number = opts$chunk_number
npeaks = length(rownames(expr))
len_chunk = round(npeaks/n_chunks)
start = (chunk_number - 1) * (len_chunk) + 1
end = start + len_chunk
if (chunk_number == n_chunks){
  end = npeaks
}

set.seed(1234)
peaks_to_run = rownames(expr)[start:end]
print(glue("selecting {n_chunks} out of {npeaks}"))
dflist = pbmclapply(peaks_to_run, mc.cores = opts$n_threads, function(i){glmnb(i, fmla, null_fmla, fmla_random, null_fmla_random)})
res = do.call(rbind, dflist)

print(res)
write.table(res, sep='\t', quote = FALSE, file = glue("{opts$prefix}.tsv"))