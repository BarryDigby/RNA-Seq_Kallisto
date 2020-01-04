#!/usr/bin/env nextflow

params.transcriptome = "Ref/Homo_sapiens.GRCh38.cdna.all.fa"
transcriptome_fasta = file( params.transcriptome )
	
params.reads = "trimmed_reads/*_r{1,2}.trimmed.fastq.gz"
Channel 
	.fromFilePairs( params.reads )
	.set { read_ch }

process Index_Transcriptome{
	publishDir "Ref/", mode:'copy'

	input:
	file transcriptome_fasta

	output:
	file "GRCh38.cDNA.idx" into indexed_transcriptome

	script:
	"""
	kallisto index -i GRCh38.cDNA.idx ${transcriptome_fasta}
	"""
}

process Quantification {
	publishDir "Kallisto_Quant/", mode:'copy'

	input:
	set val(key), file(reads) from read_ch
	file index from indexed_transcriptome

	output:
	file "*" into kallisto_out_dirs

	script:
	"""
	export BN=`cut -d'_' -f1,2,3<<<${key}`

	kallisto quant -i ${index} -t ${task.cpus} -o \$BN --bias ${reads}
	"""
}
