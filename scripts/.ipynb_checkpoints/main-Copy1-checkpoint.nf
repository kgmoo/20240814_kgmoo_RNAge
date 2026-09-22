#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'
nextflow.enable.dsl=2

// Make a channel out of an array of two strings, one a path to some set of files, one a 
def channel_glob(p) {
	chan = Channel.fromPath(p[0], type: 'any').map{it -> [it.name.replaceAll(p[1], ""), it]}
	return chan
}

// clusters = params.clusters
// cluster_colors = Channel.fromPath(params.cluster_colors, checkIfExists: true)

workflow {
	atac = channel_glob(params.atac_counts).map{it -> [it[0], it[1]]}
    //atac.view()

	//peaks = channel_glob(params.merged_peaks).map{it -> [it[0].replaceAll("pf.", ""), it[1]]}

	ch = []
	for (clust in params.clusters){
		n_chunks = params.cluster_splits[clust]
		ch << [clust, n_chunks, (1..n_chunks).toList()]
	}
	clustch = Channel.fromList(ch).transpose()
	contrasts = Channel.fromList(params.contrasts)
	
	//input = convert.out
	//.combine(peaks, by:0) // type counts peaks
	//.combine(clustch, by: 0) // type counts peaks cluster nchunks chunk
	//.combine(contrasts) // type counts peaks cluster nchunks chunk contrast
	input = atac
    .combine(clustch, by: 0)
    .combine(contrasts)
	//input.view()
	lmm(input.filter{it[4]=="age"})
    // just takes two jobs an runs them as a test
	//.take(1)

}

process lmm {
	/* NBGLMM */
	publishDir "${params.results}/lmm", mode: "rellink"
	time '48h' 
	memory '96G'
	cpus 16
	errorStrategy 'finish'
	container params.container_glmm
	tag "${cluster}.${contrast}.${n_chunks}.${chunk_number}"

	input:
	tuple val(cluster), path(sce), val(n_chunks), val(chunk_number), val(contrast)
	
	output:
	tuple val(cluster), val(contrast), val(n_chunks), path("${cluster}.${contrast}.${n_chunks}.${chunk_number}.tsv"), path("F_${cluster}.${contrast}.${n_chunks}.${chunk_number}.tsv"), path("M_${cluster}.${contrast}.${n_chunks}.${chunk_number}.tsv")

	"""
	glmnb.R --cluster $cluster --contrast $contrast --sce $sce --n_chunks $n_chunks --chunk_number $chunk_number \
	--prefix ${cluster}.${contrast}.${n_chunks}.${chunk_number} --n_threads 16
	"""
}

workflow.onComplete {
	if (workflow.success){
		subject = "ATAC glmm complete"
	}
	else {
		subject = "ATAC glmm error"
	}

	recipient = params.email

	['mail', '-s', subject, recipient].execute() << """

	Pipeline execution summary
	---------------------------
	Completed at: ${workflow.complete}
	Duration	: ${workflow.duration}
	Success		: ${workflow.success}
	workDir		: ${workflow.workDir}
	exit status : ${workflow.exitStatus}
	Error report: ${workflow.errorReport ?: '-'}
	"""
}

