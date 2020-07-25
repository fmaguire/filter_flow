#!/usr/bin/env nextflow

params.reads_pe = "$baseDir/metagenome/*{1,2}.fq"
//params.reads_comb = "$baseDir/metagenome/*comb.fq"
params.run_params = "$baseDir/run_params.csv"
params.amr_database = "$baseDir/database/nucleotide_fasta_protein_homolog_model.fasta"

// get the input read pairs
Channel
    .fromFilePairs( params.reads_pe )
    .ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
    .set{ metagenome_pair }



// parse the input csv with all the tools and params to try
Channel
    .fromPath(params.run_params)
    .splitCsv(header:true)
    .map{ row -> tuple(row.tool, row.label, row.run_params) }
    .set{ runs_ch }

// filter to just bowtie2
runs_ch
    .filter{ it[0] == "bowtie2" }
    .set{ bowtie2_params }

//// bowtie2
process prepare_bowtie2_index {
    input:
        path amr_ref from params.amr_database
    output:
        path 'amr.index*' into bowtie2_index
    script:
        """
        bowtie2-build --threads ${task.cpus} ${amr_ref} amr.index
        """
}

bowtie2_params
    .combine( metagenome_pair )
    .combine( bowtie2_index.toList() ) 
    .set{ bowtie2_run_params }

process run_bowtie2_commands {
    publishDir "bowtie2", pattern: '*.bam', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    input:
        set val(tool), val(label), val(run_params), val(read_label), path(reads), path(bowtie2_index) from bowtie2_run_params
    output:
        file '*.bam' into bowtie2_bam
    script:
        """
        bowtie2 -x bowtie2_index -U $reads $run_params -p ${task.cpus} $reads | samtools view -bS - > ${label}.bam
        """
}

