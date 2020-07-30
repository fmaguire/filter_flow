#!/usr/bin/env nextflow

/*

        Take input paired fastq reads and parse them into 
        per-tool channels.

        Secondarily, convert to single ended and fasta versions and load into
        channels for the tools that need input formatted that way.

*/
// get the input read pairs
Channel
    .fromFilePairs( params.reads_pe )
    .ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
    .set{ metagenome_fastq_pair }
metagenome_fastq_pair.into { single_end_conversion; 
                             fasta_conversion;
                             BWA_input;
                             BOWTIE2_input;
                             ARIBA_input; 
                             GROOT_input }

// generate combined single ended reads and load into channels
process build_single_ended {
    input:
        set val(label), path(reads) from single_end_conversion
    output:
        path "metagenome_se.fastq" into metagenome_fastq_se
    shell:
        """
        cat ${reads[0]} ${reads[1]} > metagenome_se.fastq
        """
}
metagenome_fastq_se.into { placeholder_input; 
                           placeholder2_input}

// generate fasta reads and load into channels
process convert_to_fasta {
    input: 
        set val(label), path(reads) from fasta_conversion
    output:
        path "metagenome_reads.fasta" into metagenome_fasta
    shell:
        """
        cat ${reads[0]} ${reads[1]} | sed -n '1~4s/^@/>/p;2~4p' > metagenome_reads.fasta
        """
}
metagenome_fasta.into{ BLASTN_input; 
                      BLASTX_input}


/* 
        Parse run parameter files for sweep settings to use in 
        homology search 

*/
Channel
    .fromPath(params.run_params_csv)
    .splitCsv(header:true)
    .map{ row -> tuple(row.tool, row.label, row.run_params) }
    .branch {
        BLASTN_params: it[0] == "blastn"
        BOWTIE2_params: it[0] == "bowtie2"
        BWA_params: it[0] == "bwa"
     }
    .set{ runs_ch }


/*

    Nucleotide Methods:
        - BLASTN
        - bowtie2
        - BWA-MEM

*/

// BLASTN
process prepare_BLASTN_database {
    conda "$baseDir/conda_envs/blast.yml"
    input:
        path amr_ref from params.amr_database
    output:
        path 'amr.blastn.db.*' into BLASTN_database
    script:
        """
        makeblastdb -dbtype nucl -in $amr_ref -out amr.blastn.db
        """
}

// filter params to just blastn and combine params with input/db to get all
// iterations of BLASTN params invocation to run correctly
runs_ch.BLASTN_params
    .combine( BLASTN_input )
    .combine( BLASTN_database.toList() )
    .set{ BLASTN_run_params }

//BLASTN defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/
process run_BLASTN_commands {
    publishDir "results/nt/blastn", pattern: '*.out6', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/blast.yml"
    input:
        set val(tool), val(label), val(run_params), path(read_fasta), path(blastn_database) from BLASTN_run_params
    output:
        file '*.out6' into BLASTN_output
    script:
        """
        blastn -query $read_fasta $run_params -num_threads ${task.cpus} -db amr.blastn.db -outfmt 6 > ${label}.out6
        """
}

// BOWTIE2
process prepare_bowtie2_index {
    conda "$baseDir/conda_envs/bowtie2.yml"
    input:
        path amr_ref from params.amr_database
    output:
        path 'amr.bowtie2.db*' into BOWTIE2_database
    script:
        """
        bowtie2-build --threads ${task.cpus} ${amr_ref} amr.bowtie2.db
        """
}

// filter to just bowtie2
runs_ch.BOWTIE2_params
    .combine( BOWTIE2_input )
    .combine( BOWTIE2_database.toList() )
    .set{ BOWTIE2_run_params }

process run_bowtie2_commands {
    publishDir "results/nt/bowtie2", pattern: '*.bam', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/bowtie2.yml"
    input:
        set val(tool), val(label), val(run_params), val(read_label), path(reads), path(bowtie2_index) from BOWTIE2_run_params
    output:
        file '*.sam' into BOWTIE2_output
    script:
        """
        bowtie2 -x amr.bowtie2.db $run_params -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} > {label}.sam
        """
        //#| samtools view -bS - > ${label}.bam
}


// BWA-MEM
process prepare_bwa_index {
    conda "$baseDir/conda_envs/bwa.yml"
    input:
        path amr_ref from params.amr_database
    output:
        path 'amr.bwa.db*' into BWA_database
    script:
        """
        bwa index -p amr.bwa.db ${amr_ref} 
        """
}

runs_ch.BWA_params
    .combine( BWA_input )
    .combine( BWA_database.toList() )
    .set{ BWA_run_params }

process run_bwa_commands {
    publishDir "results/nt/bowtie2", pattern: '*.bam', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/bwa.yml"
    input:
        set val(tool), val(label), val(run_params), val(read_label), path(reads), path(bwa_index) from BWA_run_params
    output:
        file '*.sam' into BWA_output
    script:
        """
        bwa mem -t ${task.cpus} $run_params amr.bwa.db ${reads[0]} ${reads[1]} > ${label}.sam
        """
        //#| samtools view -bS - > ${label}.bam
}



