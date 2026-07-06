
##  Functions included: LRpath, LRpath_sig, get_hs_homologs, rna_enrich_plot, and rna_enrich

LRpath<-function(file, sigvals=NULL, geneids=NULL, species, direction=NULL, min.g=10, max.g=NA, sig.cutoff=0.05,
	database="GOBP", conceptList=NULL, nullsetList=NULL,  odds.min.max=c(0.001,0.5), results_file="LRpath_results.txt") {

###########################################################################################
#  Function for LRpath: testing functional enrichment with logistic regression
#  Written by: Maureen Sartor, University of Michigan, 2009
###########################################################################################
##
##  This function uses logistic regression to test for enriched (or depleted) biological
##  categories in gene expression data. Our method models the probability 
##  of a randomly selected gene belonging to a specific category given the
##  significance level of that gene.  For categories significantly affected by
##  the experimental condition, this probability will increase as the significance
##  statistic increases.  Categories with significant p-values and positive slope 
##  coefficients are enriched with differentially expressed genes. 
##
##  Please acknowledge your use of LRpath in publications by referencing:
##  Sartor MA, Leikauf GD, and Medvedovic M.  LRpath: a logistic regression approach for
##  identifying enriched biological groups in gene expression data. Bioinformatics. 
##  25(2):211-7, 2009.
##
##  Inputs:
##  9 parameters: sigvals, geneids, species, direction, min.g, max.g, sig.cutoff, database, odds.min.max
##  "sigvals" a vector of p-values, same length and order as "geneids"
##  "geneids" a vector of Entrez gene IDs, may contain duplicates and missing values, same length and order as sigvals.
##  "species"  Used only for KEGG, GO, and Cytoband.  human="hsa", mouse="mmu", rat="rno", etc. (no default)
##		All other databases must be human, or first converted to human homologs
##		Other species for GO & cytoband are as follows: "dme" for fly, "dre" for zebrafish, "cel" for C.elegans, "sce" for yeast
##  "direction"  Optional parameter for directional testing.  If desired, a vector, of same length and order as sigvals.
##		Values should be negative if down-regulated and positive if up-regulated. (e.g. -1 or 1) Only the sign is used.
##  "min.g"  The minimum number of unique gene IDs analyzed in category to be tested
##		default = 10
##  "max.g"  The maximum number of unique gene IDs analyzed in category to be tested
##		default = NA (99999)
##  "sig.cutoff" Entrez gene IDs in each category with p-values<sig.cutoff will be returned
##		default = 0.05
##  "database" database to be tested- choices are "GO", "KEGG", "OMIM", "DrugBank", "miRBase",
##		"Panther","Biocarta", "Metabolite","MeSH","Pfam","MiMI","Cytoband","reactome", 
##		default = "GOBP"
##		"custom" for user provided gene sets
##		when "custom" is selected, conceptList is list of custom gene set IDs, 
##		nullsetList is list of genes connected to conceptList
##  "odds.min.max" lower and upper p-values to be used for odds ratio calculation.
##		default= c(0.001, 0.5)
##	"conceptList" List of concepts and genelist
##	"nullsetList" List of genes for external data source
##  "results_file" Name of file of output, default is 'LRpath_results.txt'

##
##  Outputs:
##  object is dataframe with the following columns:
##	1) GO, KEGG, or Concept ID	
##	2) GO, KEGG, or Concept term - name of category
##	3) Ontology (only for GO) - BP, MF, or CC
##	4) n.genes - Number of unique Entrez Gene IDs in category
##	5) coeff - coefficient of slope (positive values indicate enrichment or up-regulation; negative values 
##		   indicate depletion or down-regulation)
##      6) odds.ratio -  Odds ratio, as measure of strenth of enrichment (or depletion)
##      7) direction - "Enriched" or "Depleted", or if direction is used "Up-regulated" or "Down-regulated"
## 	8) p.value	-  P-value that slope does not equal zero, i.e. that term is enriched (or depleted)
##	9) FDR	- False Discovery Rate (Benjamini & Hochberg, 1995)
##  	10) sig.genes	- comma separated Entrez gene ids in category with p-value<"sig.cutoff"

##
##    Examples:
##    pvalues <- IBMT.results$IBMT.p[,1]
##    entrez.geneid <- unlist(mget(featureNames(affy.eset),env=hgu133aENTREZID))
##    GO.results <- LRpath(sigvals=pvalues, geneids=entrez.geneid, species="hsa", database="GO")
##    KEGG.results <- LRpath(sigvals=pvalues, geneids=entrez.geneid, database="KEGG",
##		species="hsa")
##
##########################################################################################
##    Load input from text file
##    input files should contain three columns: genename(or geneID), sigvalues(e.g. p-value), direction(e.g. logFC).
if(is.null(file)){
	stop("File does not exist!")
}

temp<-read.table(file,header=TRUE)
#sigvals=NULL, geneids=NULL,direction=NULL
if(ncol(temp)==2){
  sigvals=temp[,c(2)]
  geneids=temp[,c(1)]
}
if(ncol(temp)>=2){
  sigvals=temp[,c(2)]
  geneids=temp[,c(1)]
  direction=temp[,c(3)]
  direction=replace(direction,which(direction<0),-1)
  direction=replace(direction,which(direction>0),1)
}

##########################################################################################
library("annotate")
library("stats")
	
 # Get entire gene list using GO.db, for RNA-Enrich, must test each branch of GO separately for later weighting step to work correctly
		if (database %in% c('GOBP','GOCC','GOMF')) {
			 library("GO.db")
			 if (species=="hsa") {
				library("org.Hs.eg.db")
				xx<-as.list(org.Hs.egGO2ALLEGS)
				ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="mmu") {
			library("org.Mm.eg.db")
				xx<-as.list(org.Mm.egGO2ALLEGS)
				ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="rno") {
			library("org.Rn.eg.db")
			xx<-as.list(org.Rn.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="dme") {
			library("org.Dm.eg.db")
			xx<-as.list(org.Dm.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="dre") {
			library("org.Dr.eg.db")
			xx<-as.list(org.Dr.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="cel") {
			library("org.Ce.eg.db")
			xx<-as.list(org.Ce.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="sce") {
			library("org.Sc.sgd.db")
			xx<-as.list(org.Sc.sgdGO2ALLORFS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }

     ##Remove GO ids that aren't mapped to any Gene id
     xx<-xx[!is.na(xx)]

	## Determine which Entrez Gene IDs are annotated anywhere in GO
	## Gene ids not annotated in GO are left out of the environment
     if (species %in% c("sce")) {	
		nullset<-unique(unlist(xx)) 
	} else {
        nullset<-as.numeric(unique(unlist(xx)))
	}
}
#########  KEGG works for human, mouse, rat, and any other species that uses Entrez geneID
if (database=="KEGG") {	

	######################## Update on 06/01/2023 for not function KEGG library ##############
	library(limma)
	library(dplyr)
	GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = species, convert = TRUE)
	path_gene_df <- group_by(GeneID.PathID, PathwayID) %>% summarise(GeneID = list(GeneID))
	path_gene_df<-path_gene_df[length(path_gene_df$GeneID)>0,]

	xx = apply(path_gene_df, 1, function(x){
		x[["GeneID"]]
	})

	names(xx) = path_gene_df$PathwayID
	## Determine which Entrez Gene IDs are annotated anywhere in KEGG
	if (species %in% c("dme","sce")) {						
	   nullset <- unlist(GeneID.PathID$GeneID)	
	} else{
	nullset<- as.numeric(GeneID.PathID$GeneID)	
	nullset<- nullset[!is.na(nullset)]
	}


}

if (database=="External") {
	xx<-conceptList
	nullset<- nullsetList
	nullset<- nullset[!is.na(nullset)]
}

#######  For this, we'll use data in species specific R packages, e.g. org.Hs.eg.db for human
	if (database=="Cytoband") {
		 if (species=="hsa") {
		library("org.Hs.eg.db")
			###  Get GO list of Entrez genes
		xx1<-as.list(org.Hs.egMAP)
		 }
		 else if (species!="hsa") {
			stop("Only human is supported in Cytoband")
		 }
	
	## Currently, list names are Entrez IDs, and values are cytobands.  Need to switch.
	xx1<-xx1[!is.na(xx1)]
	cyto.1<-as.vector(sapply(xx1,function(x) { x[1] }))  ## names(xx1) are the entrez IDs to match this
	xx<-list()
	for (i in 1:length(unique(cyto.1))) {
		xx[[i]]<- names(xx1)[cyto.1==unique(cyto.1)[i]]
	}
	names(xx)<-as.character(unique(cyto.1))

	## Get null set vector for Cytoband
	nullset<- as.numeric(names(xx1))
	nullset<- nullset[!is.na(nullset)]
}

	if (database=="reactome") {
		if (species=="hsa") {
			load("/usr/share/lrpath/Constant/reactome_lrpath_hsa.Rdata")
		}
		else if (species=="mmu") {
			load("/usr/share/lrpath/Constant/reactome_lrpath_mmu.Rdata")		 
		}
		else if (species=="rno") {
			load("/usr/share/lrpath/Constant/reactome_lrpath_rno.Rdata")		 
		}
		else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for reactome")
			stop
		}
	}
		 
#---added 5/12/15, user inputs custom gene set list
	if (database=="custom"){
    if(!is.null(conceptList) & !is.null(nullsetList) ){
		xx = by(nullsetList, conceptList, list)
		xx = xx[!is.na(xx)]
		nullset = unique(as.numeric(unlist(xx)))
	}else{
				stop("Need conceptList and nullsetList to proceed!")
	}}
	
#----------------------------------------------------------------***
#---added 11/5/2020, chipenrich gene set list
	chipenrich_geneset_path = "C:/Users/wangdaha/Documents/chipenrichData/chipenrich.data/data/"
	chipenrich_set_used = FALSE
	if (database=="biocarta_pathway"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.biocarta_pathway.hsa.RData",sep = ""))
			xx = as.list(geneset.biocarta_pathway.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.biocarta_pathway.mmu.RData",sep = ""))
			xx = as.list(geneset.biocarta_pathway.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.biocarta_pathway.rno.RData",sep = ""))
			xx = as.list(geneset.biocarta_pathway.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for biocarta pathway")
			stop
		}
	}

	if(database=="ctd"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.ctd.hsa.RData",sep = ""))
			xx = as.list(geneset.ctd.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.ctd.mmu.RData",sep = ""))
			xx = as.list(geneset.ctd.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu')){
			message("only human and mouse are currently supported for CTD data")
			stop
		}
	}

	if(database=="drugbank"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.drug_bank.hsa.RData",sep = ""))
			xx = as.list(geneset.drug_bank.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.drug_bank.mmu.RData",sep = ""))
			xx = as.list(geneset.drug_bank.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.drug_bank.rno.RData",sep = ""))
			xx = as.list(geneset.drug_bank.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for DrugBank data")
			stop
		}
	}

	if(database=="hallmark"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.hallmark.hsa.RData",sep = ""))
			xx = as.list(geneset.hallmark.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for hallmark data")
			stop
		}
	}

	if(database=="immunologic"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.immunologic.hsa.RData",sep = ""))
			xx = as.list(geneset.immunologic.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for immunologic data")
			stop
		}
	}

	if(database=="mesh"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.mesh.hsa.RData",sep = ""))
			xx = as.list(geneset.mesh.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.mesh.mmu.RData",sep = ""))
			xx = as.list(geneset.mesh.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.mesh.rno.RData",sep = ""))
			xx = as.list(geneset.mesh.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for MeSH data")
			stop
		}
	}

	if(database=="metabolite"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.metabolite.hsa.RData",sep = ""))
			xx = as.list(geneset.metabolite.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.metabolite.mmu.RData",sep = ""))
			xx = as.list(geneset.metabolite.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.metabolite.rno.RData",sep = ""))
			xx = as.list(geneset.metabolite.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for metabolite data")
			stop
		}
	}

	if(database=="panther_pathway"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.panther_pathway.hsa.RData",sep = ""))
			xx = as.list(geneset.panther_pathway.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.panther_pathway.mmu.RData",sep = ""))
			xx = as.list(geneset.panther_pathway.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.panther_pathway.rno.RData",sep = ""))
			xx = as.list(geneset.panther_pathway.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for panther data")
			stop
		}
	}

	if(database=="pfam"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.pfam.hsa.RData",sep = ""))
			xx = as.list(geneset.pfam.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.pfam.mmu.RData",sep = ""))
			xx = as.list(geneset.pfam.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.pfam.rno.RData",sep = ""))
			xx = as.list(geneset.pfam.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for pfam data")
			stop
		}
	}

	if(database=="microrna"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.microrna.hsa.RData",sep = ""))
			xx = as.list(geneset.microrna.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for microrna data")
			stop
		}
	}

	if(database=="protein_interaction_biogrid"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.protein_interaction_biogrid.hsa.RData",sep = ""))
			xx = as.list(geneset.protein_interaction_biogrid.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for protein interaction biogrid data")
			stop
		}
	}

	if(database=="oncogenic"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.oncogenic.hsa.RData",sep = ""))
			xx = as.list(geneset.oncogenic.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for oncogenic data")
			stop
		}
	}

	if(database=="reactome_chipenrich"){
		chipenrich_set_used = TRUE
		if(species=="dme"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.dme.RData",sep = ""))
			xx = as.list(geneset.reactome.dme@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="dre"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.dre.RData",sep = ""))
			xx = as.list(geneset.reactome.dre@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.hsa.RData",sep = ""))
			xx = as.list(geneset.reactome.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.mmu.RData",sep = ""))
			xx = as.list(geneset.reactome.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.rno.RData",sep = ""))
			xx = as.list(geneset.reactome.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('dme','dre','hsa','mmu','rno')){
			message("only fruit fly, zebrafish, human, mouse, and rat are currently supported for reactome (chipenrich) data")
			stop
		}
	}

	if(database=="transcription_factors"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.transcription_factors.hsa.RData",sep = ""))
			xx = as.list(geneset.transcription_factors.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.transcription_factors.mmu.RData",sep = ""))
			xx = as.list(geneset.transcription_factors.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.transcription_factors.rno.RData",sep = ""))
			xx = as.list(geneset.transcription_factors.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for transcription_factors data")
			stop
		}
	}
#----------------------------------------------------------------***
sigvals[sigvals==0]<- 10^(-15)
geneids2<-geneids[!is.na(geneids)&!is.na(sigvals)&geneids %in% nullset]
sigvals2<-sigvals[!is.na(geneids)&!is.na(sigvals)&geneids %in% nullset]

uniqids<-unique(geneids2)
numuniq<-length(uniqids)
LOR.mult<-log(odds.min.max[2])-log(odds.min.max[1])

################  Directionality test?
if (!is.null(direction)) {
	direction2<-direction[!is.na(geneids)&!is.na(sigvals)&geneids %in% nullset]
}

################ Average -log(p-values) for duplicate locusids
nlp<-(-1)*log(sigvals2)
if (!is.null(direction)) {
	nlp[direction2<0]<-(-1)*nlp[direction2<0]
}
newp<-NA

if(length(uniqids)==length(geneids2)){
	newp = nlp[match(geneids2, uniqids)]
}else{
	newp = sapply(uniqids, function(x){ mean( nlp[which(geneids2==x)]) })
	newp = as.numeric(newp)
}

### Limit testing to those categories with at least min.g genes TOTAL
catsizes<-sapply(xx,length)
yy<-xx[catsizes>=min.g]
ncats<-length(yy)

siggenes <- uniqids[exp(-newp)<sig.cutoff]
if (!is.null(direction)) {
	siggenes<-ifelse(newp>0,uniqids[exp(-newp)<sig.cutoff],uniqids[exp(newp)<sig.cutoff])
}

######################################################################
###  Test all categories to which at least min ANALYZED genes belong

LRcoeff<-NA; LRpval<-NA; cattots<-NA; yyind<-NA; concept_ids<-NA
catsigIDs<-NA
ind<-0
if (is.na(max.g)) {
	max.g<-99999
}
for (i in (1:ncats)){
	if (species %in% c("dme","sce")) {							
	   catgenes <- as.vector(sapply(yy[[i]], as.character))	
	} else{							
	catgenes<-as.vector(sapply(yy[[i]],as.integer)) 
	}
	catgenes<-unique(catgenes)
	if (species %in% c("dme","sce")) {
	   catpoprows<-match(catgenes, as.vector(uniqids, mode="character"))
	} else{
	catpoprows<-match(catgenes,as.vector(uniqids,mode="integer"))
	}
	catpoprows<-catpoprows[!is.na(catpoprows)]
	

	cattot<-length(catpoprows)	#Number of analyzed genes in category
	if (cattot>=min.g&cattot<=max.g) {
		ind<-ind+1
		cat<-rep(0,numuniq)
		cat[catpoprows]<- 1
		forLR<- as.data.frame(cbind(cat,newp))
		names(forLR)<-c("cat","nlogp")

		glm.lrGO <- glm(cat ~ nlogp, family=binomial(link="logit"),forLR)
		lrGO<-summary(glm.lrGO)

		###########  Extract p-value and coefficient
		LRcoeff[ind]<-lrGO$coefficients["nlogp","Estimate"]
		LRpval[ind]<-lrGO$coefficients["nlogp","Pr(>|z|)"]
		cattots[ind]<-cattot
		#yyind[ind]<-i
		concept_ids[ind]<-names(yy[i])

		catsig.IDlist<-intersect(siggenes,catgenes)
		if (!is.null(direction)) {
		   if (LRcoeff[ind]>0) {
			catsig.IDlist<-intersect(catsig.IDlist,uniqids[newp>0])
		   }
		   else { catsig.IDlist<-intersect(catsig.IDlist,uniqids[newp<0]) }
		}
		catsigIDs[ind]<-paste(catsig.IDlist[order(catsig.IDlist)],collapse=", ")
		catsigIDs[ind]<- ifelse(length(catsig.IDlist)==0,NA,catsigIDs[ind]) #added 7/08
      }
	  # print(i)
	if ((i-ncats/10)>0&(i-ncats/10)<1.1) {
		print("10% categories finished.")
	}
	if ((i-ncats/5)>0&(i-ncats/5)<1.1) {
		print("20% categories finished.")
	}
	if ((i-ncats/2)>0&(i-ncats/2)<1.1) {
		print("50% categories finished.")
	}
	if ((i-9*ncats/10)>0&(i-9*ncats/10)<1.1) {
		print("90% categories finished.")
	}

}

odds.ratio<-exp(LOR.mult*LRcoeff)
	#####calculate FDR
	BenjFDR<-p.adjust(LRpval,"BH")
		
if (database %in% c("GOBP","GOCC","GOMF")) {
	############## Annotate with GO names and Ontology
	GOterms<-as.list(GOTERM)
	GOterms<-GOterms[!is.na(GOterms)]
	terms<-as.vector(sapply(GOterms,Term))  
	Ontologies<-as.vector(sapply(GOterms,Ontology))
	AllGOids<-as.vector(names(GOterms))
	terms2<-cbind(AllGOids,terms,Ontologies)

	GOrows<-match(concept_ids,AllGOids)
	GOannot<-terms2[GOrows,]
	LRresults<-cbind(concept_ids,GOannot[,c("terms","Ontologies")],cattots,LRcoeff,
	     odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	names(LRresults)<- c("Concept.ID", "Concept.name", "Ontology", "n.genes", "coeff", "odds.ratio",
		"p.value", "FDR", "sig.genes")
}

if (database=="KEGG") {
	# KEGGterms<- as.list(KEGGPATHNAME2ID)
	# kterms<- as.vector(names(KEGGterms))
	# k.ids<- paste(species,unlist(KEGGterms),sep="")
	# keggrows<- match(concept_ids,k.ids)
	# kegg.annot<-kterms[keggrows]
	# LRresults<-cbind(concept_ids,kegg.annot,cattots,LRcoeff,
	#      odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	# names(LRresults)<- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
	# 	"p.value", "FDR", "sig.genes")

	######################## Update on 06/01/2023 for not function KEGG library ##############
	  KEGGterms <- limma::getKEGGPathwayNames(species.KEGG = species, 
                                         remove.qualifier = TRUE)
	  kterms <- as.vector(KEGGterms$Description)
	  k_ids <- paste("path:", KEGGterms$PathwayID, sep = "")
	  keggrows <- match(concept_ids, k_ids)
	  kegg.annot <- kterms[keggrows]

	  LRresults <- cbind(concept_ids, kegg.annot, cattots, LRcoeff,
	                     odds.ratio, as.data.frame(LRpval), BenjFDR, catsigIDs)
	  names(LRresults) <- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
	                        "p.value", "FDR", "sig.genes")
}
if (database=="Cytoband") {
	#Cytobandrows = match(concept_ids, Cytoband.annot[, "CONCEPT_ID"])	
	LRresults<-cbind(concept_ids,concept_ids,cattots,LRcoeff,
	     odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	names(LRresults)<- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
		"p.value", "FDR", "sig.genes")
}

if (database=="External") {
	LRresults<-cbind(concept_ids,concept_ids,cattots,LRcoeff,
	     odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	names(LRresults)<- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
		"p.value", "FDR", "sig.genes")
}

if (database=="reactome") {
	reactome.annot = reactome_desc$reactome_name[match(concept_ids, reactome_desc$reactome_id)]
	LRresults<-cbind(concept_ids,reactome.annot,cattots,LRcoeff,
	     odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	names(LRresults)<- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
		"p.value", "FDR", "sig.genes")
}

if (database=="custom") {
	LRresults<-cbind(concept_ids,annot=concept_ids,cattots,LRcoeff,
	     odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	names(LRresults)<- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
		"p.value", "FDR", "sig.genes")
}

if(chipenrich_set_used == TRUE){
	LRresults<-cbind(concept_ids,concept_ids,cattots,LRcoeff,
	     odds.ratio,as.data.frame(LRpval),BenjFDR,catsigIDs)
	names(LRresults)<- c("Concept.ID", "Concept.name", "n.genes", "coeff", "odds.ratio",
		"p.value", "FDR", "sig.genes")
}

LRresults=LRresults[!is.na(LRresults$p.value),];
LRresults = LRresults[order(LRresults$FDR),]

write.table(LRresults,file=results_file,sep="\t",row.names=F,quote=F)
LRresults

} #End function LRpath


LRpath_sig<-function(LRresults, n = 10) {
	LRsiggenes = LRresults[order(LRresults[, "p.value"])[1:n],]

LRsiggenes
}
	

get_hs_homologs <- function(geneids) {
	library(homolog.db)
	human.ids<- as.list(homologTAXID)
	human.ids<-names(human.ids)[human.ids==9606] 
	print("Human IDs acquired.")
	xx<-as.list(homologHOMOLOG)	
	hids<- xx[match(geneids,names(xx))]
	homologs <- names(xx)[names(xx) %in% human.ids][match(hids,xx)]

homologs
}


rna_enrich_plot = function(sigvals,geneids,avg_readcount,filename="rna_enrich_plot.jpg", plot_height=480, plot_width=480){
	library(mgcv); library(rms);
	de_df = data.frame(geneid = geneids, p.value = sigvals, readcount = avg_readcount)
	de_df = na.omit(de_df)
	de_df = de_df[de_df$readcount>=5,]
	de_df = de_df[order(de_df$p.value),]
		
	rownames(de_df) = 1:nrow(de_df)
	de_df$group = ceiling(as.numeric(rownames(de_df))/25)
	gavgrc = as.numeric(by(de_df$readcount, de_df$group, mean))
	gavgp = as.numeric(by(de_df$p.value, de_df$group, mean))
	log10_gavgrc = log10(gavgrc)
	log10_gavgp = -log10(gavgp + min(gavgp[gavgp!=0]))

	de_df$log10_p.value = (-1)*log10(de_df$p.value + min(de_df$p.value[de_df$p.value !=0],na.rm=T))
	de_df$log10_readcount = log10(de_df$readcount)
	de_df = de_df[!is.na(de_df$log10_p.value),]
	model = "log10_p.value ~ s(log10_readcount,bs='cr')"
	  # Compute binomial spline fit.  
	fit = gam(as.formula(model),data=de_df,family="gaussian")
	  # Compute weights for each gene, based on the predicted prob(peak) for each gene. 
	pP = fitted(fit);
	 
	jpeg(filename, height=plot_height, width=plot_width)
	par(mar= c(5, 5, 4, 2) + 0.1)
	  
	plot(log10_gavgp ~ log10_gavgrc, ylab='-log10(p-value) - binned 25 genes', xlab='log10(avg read count) - binned 25 genes',pch=20, cex.axis=1.5, cex.lab=1.5) 
	lines(pP[order(de_df$log10_readcount)]~sort(de_df$log10_readcount),col='red',lwd=2)
	
	dev.off()
}

#RNA-Enrich test	

rna_enrich = function(file, sigvals=NULL, geneids=NULL, avg_readcount, species, direction=NULL, min.g = 10,  max.g=NA, sig.cutoff=0.05, database="GOBP", adj_readcount=TRUE, read_lim=5,conceptList=NULL, nullsetList=NULL, filename="rna_enrich_plot.jpg",plot_height=480, plot_width=480,results_file="RNA-Enrich_results.txt") {
	#input file should have at least three columns, i.e. geneid, p.value, readcount. 
	#The fourth column can be the direction, i.e. log2foldChange

	if(is.null(file)){
	stop("File does not exist!")
	}

	df = read.table(file,header=TRUE)
	
	if(ncol(df) == 3){
		geneids = df[,1]
		sigvals = df[,2]
		avg_readcount = df[,3]
		de_df = data.frame(geneid = geneids, p.value = sigvals, readcount = avg_readcount)
		de_df = stats::aggregate(cbind(p.value,readcount) ~ geneid, de_df, mean)
	}
	if(ncol(df) == 4){
		geneids = df[,1]
		sigvals = df[,2]
		avg_readcount = df[,3]
		direction = df[,4]
		de_df = data.frame(geneid = geneids, p.value = sigvals, readcount = avg_readcount)
		de_df$direction = direction
		de_df = stats::aggregate(cbind(p.value,readcount,direction) ~ geneid, de_df, mean) 
		de_df$direction = ifelse(de_df$direction>0,1,-1)
	}
	# de_df = data.frame(geneid = geneids, p.value = sigvals, readcount = avg_readcount)

	# if(!is.null(direction)){
	# 	de_df$direction = direction
	# 	de_df = stats::aggregate(cbind(p.value,readcount,direction) ~ geneid, de_df, mean) 
	# 	de_df$direction = ifelse(de_df$direction>0,1,-1)
	# }else{ 	
	# 	de_df = stats::aggregate(cbind(p.value,readcount) ~ geneid, de_df, mean)
	# }
		
	rna_enrich_plot(sigvals=sigvals, geneids=geneids, avg_readcount = avg_readcount, filename=filename, plot_height=plot_height, plot_width=plot_width)

	
	de_df = na.omit(de_df)
	de_df = subset(de_df, geneid !='')
	
	if(!species %in% c('dme','sce')){
		de_df$geneid = suppressWarnings(as.numeric(as.character(de_df$geneid)))
		de_df = na.omit(de_df)
	}

	
library("annotate")
library("stats")
	
 # Get entire gene list using GO.db, for RNA-Enrich, must test each branch of GO separately for later weighting step to work correctly
		if (database %in% c('GOBP','GOCC','GOMF')) {
			 library("GO.db")
			 if (species=="hsa") {
				library("org.Hs.eg.db")
				xx<-as.list(org.Hs.egGO2ALLEGS)
				ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="mmu") {
			library("org.Mm.eg.db")
				xx<-as.list(org.Mm.egGO2ALLEGS)
				ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="rno") {
			library("org.Rn.eg.db")
			xx<-as.list(org.Rn.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="dme") {
			library("org.Dm.eg.db")
			xx<-as.list(org.Dm.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="dre") {
			library("org.Dr.eg.db")
			xx<-as.list(org.Dr.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="cel") {
			library("org.Ce.eg.db")
			xx<-as.list(org.Ce.egGO2ALLEGS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }
			 else if (species=="sce") {
			library("org.Sc.sgd.db")
			xx<-as.list(org.Sc.sgdGO2ALLORFS)
			ontols = Ontology(names(xx))
				if(database == 'GOBP'){
					xx = xx[ontols=='BP']
				}
				if(database == 'GOCC'){
					xx = xx[ontols=='CC']
				}
				if(database == 'GOMF'){
					xx = xx[ontols=='MF']
				}
			 }

     ##Remove GO ids that aren't mapped to any Gene id
     xx<-xx[!is.na(xx)]

	## Determine which Entrez Gene IDs are annotated anywhere in GO
	## Gene ids not annotated in GO are left out of the environment
     if (species %in% c("sce")) {	
		nullset<-unique(unlist(xx)) 
	} else {
        nullset<-as.numeric(unique(unlist(xx)))
	}
}
#########  KEGG works for human, mouse, rat, and any other species that uses Entrez geneID
if (database=="KEGG") {	
	######################## Update on 06/01/2023 for not function KEGG library ##############
	library(limma)
	library(dplyr)
	GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = species, convert = TRUE)
	path_gene_df <- group_by(GeneID.PathID, PathwayID) %>% summarise(GeneID = list(GeneID))
	path_gene_df<-path_gene_df[length(path_gene_df$GeneID)>0,]

	xx = apply(path_gene_df, 1, function(x){
		x[["GeneID"]]
	})

	names(xx) = path_gene_df$PathwayID
	## Determine which Entrez Gene IDs are annotated anywhere in KEGG
	if (species %in% c("dme","sce")) {						
	   nullset <- unlist(GeneID.PathID$GeneID)	
	} else{
	nullset<- as.numeric(GeneID.PathID$GeneID)	
	nullset<- nullset[!is.na(nullset)]
	}

}

if (database=="External") {
	xx<-conceptList
	nullset<- nullsetList
	nullset<- nullset[!is.na(nullset)]
}

#######  For this, we'll use data in species specific R packages, e.g. org.Hs.eg.db for human
	if (database=="Cytoband") {
		 if (species=="hsa") {
		library("org.Hs.eg.db")
			###  Get GO list of Entrez genes
		xx1<-as.list(org.Hs.egMAP)
		 }
		 else if (species!="hsa") {
			stop("Only human is supported in Cytoband")
		 }
	
	## Currently, list names are Entrez IDs, and values are cytobands.  Need to switch.
	xx1<-xx1[!is.na(xx1)]
	cyto.1<-as.vector(sapply(xx1,function(x) { x[1] }))  ## names(xx1) are the entrez IDs to match this
	xx<-list()
	for (i in 1:length(unique(cyto.1))) {
		xx[[i]]<- names(xx1)[cyto.1==unique(cyto.1)[i]]
	}
	names(xx)<-as.character(unique(cyto.1))

	## Get null set vector for Cytoband
	nullset<- as.numeric(names(xx1))
	nullset<- nullset[!is.na(nullset)]
}

	if (database=="reactome") {
		if (species=="hsa") {
			load("/usr/share/LrpathFileStore/reactome_lrpath_hsa.Rdata")
		}
		else if (species=="mmu") {
			load("/usr/share/LrpathFileStore/reactome_lrpath_mmu.Rdata")		 
		}
		else if (species=="rno") {
			load("/usr/share/LrpathFileStore/reactome_lrpath_rno.Rdata")		 
		}
		else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for reactome")
			stop
		}
	}

#---added 5/12/15, user inputs custom gene set list
	if (database=="custom"){
    if(!is.null(conceptList) & !is.null(nullsetList) ){
		xx = by(nullsetList, conceptList, list)
		xx = xx[!is.na(xx)]
		nullset = unique(as.numeric(unlist(xx)))
	}else{
				stop("Need conceptList and nullsetList to proceed!")
	}}
	

#---added 11/11/2020
#need change this path to your chipenrich data folder
	chipenrich_geneset_path = "C:/Users/wangdaha/Documents/chipenrichData/chipenrich.data/data/"
	chipenrich_set_used = FALSE
	if (database=="biocarta_pathway"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.biocarta_pathway.hsa.RData",sep = ""))
			xx = as.list(geneset.biocarta_pathway.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.biocarta_pathway.mmu.RData",sep = ""))
			xx = as.list(geneset.biocarta_pathway.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.biocarta_pathway.rno.RData",sep = ""))
			xx = as.list(geneset.biocarta_pathway.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for biocarta pathway")
			stop
		}
	}

	if(database=="ctd"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.ctd.hsa.RData",sep = ""))
			xx = as.list(geneset.ctd.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.ctd.mmu.RData",sep = ""))
			xx = as.list(geneset.ctd.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu')){
			message("only human and mouse are currently supported for CTD data")
			stop
		}
	}

	if(database=="drugbank"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.drug_bank.hsa.RData",sep = ""))
			xx = as.list(geneset.drug_bank.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.drug_bank.mmu.RData",sep = ""))
			xx = as.list(geneset.drug_bank.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.drug_bank.rno.RData",sep = ""))
			xx = as.list(geneset.drug_bank.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for DrugBank data")
			stop
		}
	}

	if(database=="hallmark"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.hallmark.hsa.RData",sep = ""))
			xx = as.list(geneset.hallmark.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for hallmark data")
			stop
		}
	}

	if(database=="immunologic"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.immunologic.hsa.RData",sep = ""))
			xx = as.list(geneset.immunologic.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for immunologic data")
			stop
		}
	}

	if(database=="mesh"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.mesh.hsa.RData",sep = ""))
			xx = as.list(geneset.mesh.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.mesh.mmu.RData",sep = ""))
			xx = as.list(geneset.mesh.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.mesh.rno.RData",sep = ""))
			xx = as.list(geneset.mesh.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for MeSH data")
			stop
		}
	}

	if(database=="metabolite"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.metabolite.hsa.RData",sep = ""))
			xx = as.list(geneset.metabolite.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.metabolite.mmu.RData",sep = ""))
			xx = as.list(geneset.metabolite.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.metabolite.rno.RData",sep = ""))
			xx = as.list(geneset.metabolite.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for metabolite data")
			stop
		}
	}

	if(database=="panther_pathway"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.panther_pathway.hsa.RData",sep = ""))
			xx = as.list(geneset.panther_pathway.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.panther_pathway.mmu.RData",sep = ""))
			xx = as.list(geneset.panther_pathway.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.panther_pathway.rno.RData",sep = ""))
			xx = as.list(geneset.panther_pathway.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for panther data")
			stop
		}
	}

	if(database=="pfam"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.pfam.hsa.RData",sep = ""))
			xx = as.list(geneset.pfam.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.pfam.mmu.RData",sep = ""))
			xx = as.list(geneset.pfam.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.pfam.rno.RData",sep = ""))
			xx = as.list(geneset.pfam.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for pfam data")
			stop
		}
	}

	if(database=="microrna"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.microrna.hsa.RData",sep = ""))
			xx = as.list(geneset.microrna.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for microrna data")
			stop
		}
	}

	if(database=="protein_interaction_biogrid"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.protein_interaction_biogrid.hsa.RData",sep = ""))
			xx = as.list(geneset.protein_interaction_biogrid.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for protein interaction biogrid data")
			stop
		}
	}

	if(database=="oncogenic"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.oncogenic.hsa.RData",sep = ""))
			xx = as.list(geneset.oncogenic.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa')){
			message("only human is currently supported for oncogenic data")
			stop
		}
	}

	if(database=="reactome_chipenrich"){
		chipenrich_set_used = TRUE
		if(species=="dme"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.dme.RData",sep = ""))
			xx = as.list(geneset.reactome.dme@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="dre"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.dre.RData",sep = ""))
			xx = as.list(geneset.reactome.dre@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.hsa.RData",sep = ""))
			xx = as.list(geneset.reactome.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.mmu.RData",sep = ""))
			xx = as.list(geneset.reactome.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.reactome.rno.RData",sep = ""))
			xx = as.list(geneset.reactome.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('dme','dre','hsa','mmu','rno')){
			message("only fruit fly, zebrafish, human, mouse, and rat are currently supported for reactome (chipenrich) data")
			stop
		}
	}

	if(database=="transcription_factors"){
		chipenrich_set_used = TRUE
		if(species=="hsa"){
			load(paste(chipenrich_geneset_path,"geneset.transcription_factors.hsa.RData",sep = ""))
			xx = as.list(geneset.transcription_factors.hsa@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="mmu"){
			load(paste(chipenrich_geneset_path,"geneset.transcription_factors.mmu.RData",sep = ""))
			xx = as.list(geneset.transcription_factors.mmu@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if(species=="rno"){
			load(paste(chipenrich_geneset_path,"geneset.transcription_factors.rno.RData",sep = ""))
			xx = as.list(geneset.transcription_factors.rno@set.gene)
			nullset = unique(as.numeric(unlist(xx)))
			nullset = nullset[!is.na(nullset)]
		}else if (!species %in% c('hsa','mmu','rno')){
			message("only human, mouse, and rat are currently supported for transcription_factors data")
			stop
		}
	}
	#-------------------
	catsizes<-sapply(xx,length);
	if(is.na(max.g)){max.g = max(catsizes)}
	xx<-xx[catsizes>=min.g & catsizes<=max.g];
	ncats<-length(xx);

	#make empty d dataframe to initiate gpw creation
	d_allgenes = data.frame(geneid=nullset[which(!nullset %in% de_df$geneid)],p.value=NA,readcount=0);
	if(!is.null(direction)){
		d_allgenes$direction = NA;
		d_allgenes$direction = as.numeric(d_allgenes$direction)
	};
	gpw = rbind(de_df,d_allgenes);


	# Add log10 p-values to gpw, add min p[!=0] to all pvals so can log transform later
	gpw$log10_p.value = (-1)*log10(gpw$p.value + min(gpw$p.value[gpw$p.value !=0],na.rm=T));

	# remove genes with read lim readcounts (default is 5) and log readcount
	if(read_lim==0){ 
		gpw = subset(gpw, readcount>read_lim);
	}else{
		gpw = subset(gpw, readcount>=read_lim);
	}
	gpw$log10_readcount = log10(gpw$readcount);

	# Sort by readcount 
	gpw = gpw[order(gpw$readcount,decreasing=T),];

	gpw = gpw[!is.na(gpw$log10_p.value),];
	# Create model - adjust for readcount
	if(adj_readcount==TRUE){
		model = "log10_p.value ~ s(log10_readcount,bs='cr')";
	}else{model = "log10_p.value ~ 1"};

	# Compute binomial spline fit.
	library(mgcv); library(rms);
	fit = gam(as.formula(model),data=gpw,family="gaussian");

	# Compute weights for each gene, based on the predicted prob(DE) for each gene. 
	pP = fitted(fit);
	w0 = 1 / (pP/mean(gpw$log10_p.value,na.rm=T));
	w0 = w0 / mean(w0,na.rm=T);

	gpw$weight = w0;
	gpw$prob_P = pP;
	gpw$resid.dev = resid(fit,type="deviance");

	# change log pvals if directional test
	if (!is.null(direction)) {
		gpw$log10_p.value = gpw$log10_p.value * gpw$direction
	}

	cols = c("geneid","readcount","log10_readcount","p.value","log10_p.value","weight","prob_P","resid.dev");
	gpw = subset(gpw,select=cols);

	#now run test	

	lrm.fast = function(x,y) {
		fit = lrm.fit(x,y);
		vv = diag(fit$var);
		cof = fit$coef;
		z = cof/sqrt(vv);
		pval = pchisq(z^2,1,lower.tail=F);
		c(cof[2],pval[2]);
	};
	
	# Restrict our genes/weights/peaks to only those genes in the genesets. 
	gpw = subset(gpw,geneid %in% nullset);
	
	# Re-normalize weights so that mean is 1 (if any genes were removed from previous gpw step)
	gpw$weight = gpw$weight / mean(gpw$weight);
	
	r_pvals = c();
	r_effects = c();
	r_go_ids = c();
	r_go_genes = c();
	r_go_genes_num = c();
	r_go_genes_DE = c();
	r_go_genes_DE_num = c();
	r_go_genes_avg_length = c();
	for (i in 1:length(xx)) {
		
		go_id = names(xx)[i];
		
		if(species %in% c('dme','sce')){
			go_genes = as.character(xx[[i]]);
		}else{
			go_genes = as.numeric(as.character(xx[[i]]));
		}
		# Eliminate GO genes that aren't in the gpw. 
		#go_genes = Filter(function(x) x %in% gpw$geneid,go_genes);
		go_genes = go_genes[go_genes %in% gpw$geneid];
		
		# A boolean vector of gene membership of all geneids in go_genes
		b_genes = gpw$geneid %in% go_genes;
		sg_go = gpw$log10_p.value[b_genes];
		wg_go = gpw$weight[b_genes];
		
		r_go_genes[i] = paste(go_genes,collapse=";");
		r_go_genes_num[i] = length(go_genes);
		
		# go_genes_DE = gpw$geneid[b_genes][sg_go==1];
		r_go_genes_DE_num[i] = length(gpw$geneid[gpw$p.value<=sig.cutoff & b_genes==TRUE]);
		go_genes_DE = gpw$geneid[gpw$p.value<=sig.cutoff & b_genes==TRUE]; 
		if(length(go_genes_DE)==0){go_genes_DE=NA}
		r_go_genes_DE[i] = paste(go_genes_DE,collapse=", ");
		
    	
    	# The model is still essentially InGO ~ s(log10_length,bs='cr') + -log10_p.value
    	# except we're using the weights that are calculated once...
    	# Also, y must be binary, so as.numeric(b_genes) (goterm) must be 
    	# predicted against the weights...
		testm = cbind(y=as.numeric(b_genes), x=gpw$weight*gpw$log10_p.value);
		
		if(sum(testm[,'y']==1) == nrow(gpw) | sum(testm[,'y']) ==0){
		ep = c(NA,NA)
		}else{ep = lrm.fast(testm[,"x"], testm[,"y"])}
		
		r_effects[i] = ep[1];
		r_pvals[i] = ep[2];
		r_go_ids[i] = go_id;
		
		if ((i-ncats/10)>0&(i-ncats/10)<1.1) {
		print("10% categories finished.")
		};
		if ((i-ncats/5)>0&(i-ncats/5)<1.1) {
			print("20% categories finished.")
		};
		if ((i-ncats/2)>0&(i-ncats/2)<1.1) {
			print("50% categories finished.")
		};
		if ((i-9*ncats/10)>0&(i-9*ncats/10)<1.1) {
			print("90% categories finished.")
		};
	}
	
	fdr = p.adjust(r_pvals,method="BH");
  
	is_depleted = r_effects < 0; 
	is_enriched = r_effects > 0;
  
	enr = rep(NA,length(r_go_ids));
	if(!is.null(direction)){
		enr[is_depleted] = "down";
		enr[is_enriched] = "up";
	}else{
		enr[is_depleted] = "depleted";
		enr[is_enriched] = "enriched";
    }
 ##---annotate concept_id/r_go_ids with description, right now just have GO and KEGG, NAs for others
  concept_descriptions = NA
  if(database %in% c('GOBP','GOCC','GOMF')){
	concept_descriptions = Term(r_go_ids);
  };
  if(database=="KEGG") {
	# KEGGterms<- as.list(KEGGPATHNAME2ID);
	# kterms<- as.vector(names(KEGGterms));
	# k.ids<- paste(species,unlist(KEGGterms),sep="");
	# keggrows<- match(r_go_ids,k.ids);
	# concept_descriptions<-kterms[keggrows];


	######################## Update on 06/01/2023 for not function KEGG library ##############
	  KEGGterms <- limma::getKEGGPathwayNames(species.KEGG = species, 
                                         remove.qualifier = TRUE)
	  kterms <- as.vector(KEGGterms$Description)
	  k_ids <- paste("path:", KEGGterms$PathwayID, sep = "")
	  keggrows <- match(concept_ids, k_ids)
	  
	concept_descriptions<-kterms[keggrows]
  };
  if (database=="reactome") {
	concept_descriptions = reactome_desc$reactome_name[match(r_go_ids, reactome_desc$reactome_id)]
  };
  if (database=="custom") {
	concept_descriptions = r_go_ids
  };

	results = data.frame(
		"Concept ID"=r_go_ids,
		"Concept name" = concept_descriptions,
		"ontology" = database, 
		"n.genes"=r_go_genes_num, #changed 7/8, incorrect before
		"coeff"=r_effects,
		"odds.ratio"=exp(r_effects),
		"status"=enr,
		"p.value"=r_pvals,
		"FDR"=fdr,
		"sig.genes"=r_go_genes_DE,
		stringsAsFactors=F
	);
	results = results[ !is.na(results$p.value),];

	results = results[order(results$FDR,results$p.value),];
	write.table(results, file=results_file, sep="\t",row.names=F, quote=F)
	return(results);
}