#!/usr/bin/env nextflow

params.transcriptome = "/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Reference/Homo_sapiens.GRCh38.cdna.all.fa"
transcriptome_fasta = file( params.transcriptome )

read_ch = Channel.fromPath("/data/MA5112/Practicals/RNA-Seq/Kallisto_Practical/Data/*.fastq.gz")

process Index_Transcriptome{
        publishDir "Reference/", mode:'copy'

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
        file reads from read_ch
        file index from indexed_transcriptome

        output:
        file "*" into kallisto_out_dirs

        script:
        """
        export BN=`cut -d'.' -f1<<<${reads}`

        kallisto quant --single -l 200 -s 30 -i ${index} -t ${task.cpus} -o \$BN --bias ${reads}
        """
}
