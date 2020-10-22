#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*

        Take input paired fastq reads and parse them into 
        per-tool channels.
        Secondarily, convert to single ended and fasta versions and load into
        channels for the tools that need input formatted that way.

*/

include { build_single_ended; 
          convert_to_fasta; 
          convert_to_amino_acid; } from './modules/preprocessing.nf'

/*
Nucleotide Methods vs Nucleotide database:
        - BLASTN
        - BOWTIE2
        - BWA
        - groot
        - ariba
        - KMA
        - HMMSEARCHNT
*/

include { prepare_BLASTN_database; run_BLASTN_commands; 
          prepare_BOWTIE2_index; run_BOWTIE2_commands; 
          prepare_BWA_index; run_BWA_commands; 
          prepare_GROOT_database; run_GROOT_commands; 
          prepare_ARIBA_database; run_ARIBA_commands;
          prepare_KMA_database; run_KMA_commands;
          prepare_HMMSEARCHNT_database; run_HMMSEARCHNT_commands } from './modules/nt.nf'


/*
Nucleotide Methods vs Protein database:
        - BLASTX
        - DIAMONDBLASTX
        - PALADIN
        - MMSEQSBLASTX
*/

include { prepare_BLASTX_BLASTP_database; run_BLASTX_commands; 
          prepare_DIAMOND_databases; run_DIAMONDBLASTX_commands; 
          prepare_PALADIN_database; run_PALADIN_commands; 
          prepare_MMSEQS_database; run_MMSEQSBLASTX_commands} from "./modules/nt_to_aa.nf"

/*
Protein Methods vs Protein database:
        - BLASTP
        - DIAMONDBLASTP
        - MMSEQS
*/

//include { } from "./modules/aa.nf"

// * DIAMOND-BLASTP
//// database is already generated from previous diamond invocation
//runs_ch.DIAMONDBLASTP_params
//    .map{ it -> it[0,1,2] }
//    .combine( DIAMONDBLASTP_input )
//    .combine( DIAMONDBLASTP_database.toList() )
//    .set{ DIAMONDBLASTP_run_params }
//
//process run_DIAMONDBLASTP_commands {
//    tag { "DIAMONDBLASTP: ${label}" }
//    publishDir "results/aa/diamond-blastp", pattern: "${label}.out6", mode: "copy"
//    conda "$baseDir/conda_envs/diamond.yml"
//    input:
//        tuple val(tool), val(label), val(run_params), path(read_aa_fasta), path(diamond_database) from DIAMONDBLASTP_run_params
//    output:
//        path "${label}.out6" into DIAMONDBLASTP_output
//    script:
//        """
//        diamond blastp -p ${task.cpus} --max-target-seqs 1 --db ${diamond_database} --outfmt 6 --out ${label}.out6 ${run_params} -q ${read_aa_fasta}
//        """
//} 
//
//
//// * BLASTP
//runs_ch.BLASTP_params
//    .map{ it -> it[0,1,2] }
//    .combine( BLASTP_input )
//    .combine( BLASTP_database.toList() )
//    .set{ BLASTP_run_params }
//
////BLASTP defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastp_application_options/
//process run_BLASTP_commands {
//    tag { "BLASTP: ${label}" }
//    publishDir "results/aa/blastp", pattern: "${label}.out6", mode: "copy"
//    conda "$baseDir/conda_envs/blast.yml"
//    input:
//        tuple val(tool), val(label), val(run_params), path(read_aa_fasta), path(blast_database) from BLASTP_run_params
//    output:
//        path "${label}.out6" into BLASTP_output
//    script:
//        """
//        blastp -max_target_seqs 1  -query $read_aa_fasta $run_params -num_threads ${task.cpus} -db amr.prot.db -outfmt 6 > ${label}.out6
//        """
//}
//
//runs_ch.MMSEQS_blastp_params
//    .map{ it -> it[0,1,2] }
//    .combine( MMSEQS_blastp_input )
//    .combine( MMSEQS_blastp_database.toList() )
//    .set{ MMSEQS_blastsp_run_params }
//
//
//process run_MMSEQS_blastp_commands {
//    tag { "MMSEQS: ${label}" }
//    publishDir: "results/aa/mmseqs_blastp", pattern: "", mode: "copy"
//    conda: "$baseDir/conda_envs/mmseqs2.yml"
//    input:
//        tuple val(tool), val(label), val(run_params), path(read_aa_fasta), path(mmseqs_database), path(mmseqs_tmp) from MMSEQS_blastp_run_params
//    output: 
//        asd
//    script:
//        """
//        mmseqs easy-search --threads ${task.cpus} ${read_aa_fasta} ${mmseqs_database}
//        """
//}
// As all sequences are fragemtns 
//hmmsearch --fragthresh 1



workflow {
    /*
        Set up the input data for the workflow

    */

    // get the input read pairs
    fastq_paired_end_reads = Channel.fromFilePairs( params.reads_pe )
                                     .ifEmpty{ error "Cannot find any reads matching: ${params.reads_pe}" }

    // generate combined single ended reads and load into channels
    fastq_single_ended_reads = build_single_ended( fastq_paired_end_reads );

    // generate fasta reads and load into channels
    fasta_single_ended_reads = convert_to_fasta( fastq_paired_end_reads );

    // convert input metagenome to amino acid orfs using ORFM
    fasta_single_ended_amino_acid = convert_to_amino_acid( fasta_single_ended_reads );

    /* 
        Parse run parameter files for sweep settings to use in 
        homology search 

    */
    runs_ch = Channel
                .fromPath(params.run_params_csv)
                .splitCsv(header:true)
                .map{ row -> tuple(row.tool, row.label, row.run_params, row.db_params) }
                .branch {
                    BLASTN_params: it[0] == "blastn"
                    BLASTX_params: it[0] == "blastx"
                    BLASTP_params: it[0] == "blastp"
                    BOWTIE2_params: it[0] == "bowtie2"
                    BWA_params: it[0] == "bwa"
                    GROOT_params: it[0] == "groot"
    	        	ARIBA_params: it[0] == "ariba"
                    DIAMONDBLASTX_params: it[0] == 'diamond-blastx'
                    DIAMONDBLASTP_params: it[0] == 'diamond-blastp'
                    KMA_params: it[0] == 'kma'
                    HMMSEARCHNT_params: it[0] == 'hmmsearch_nt'
                    PALADIN_params: it[0] == 'paladin'
                    MMSEQSBLASTX_params: it[0] == 'mmseqsblastx'
                 }


    /*

        Nucleotide Methods vs Nucleotide database:

    */

    // BLASTN
    BLASTN_database = prepare_BLASTN_database( params.amr_nucl_database );

    // filter params to just blastn and combine params with input/db to get all
    // iterations of BLASTN params invocation to run correctly
    BLASTN_run_params = runs_ch.BLASTN_params
                                .map{ it -> it[0,1,2] }
                                .combine( fasta_single_ended_reads )
                                .combine( BLASTN_database.toList() )

    //BLASTN defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/
    BLASTN_output = run_BLASTN_commands( BLASTN_run_params );
    
    // parse blast on reads (ORFs need separate)
    //parsed_BLASTN_output = parse_BLAST_output( BLASTN_output );

    //// BOWTIE2
    BOWTIE2_database = prepare_BOWTIE2_index( params.amr_nucl_database );

    //// filter to just bowtie2
    BOWTIE2_run_params = runs_ch.BOWTIE2_params
                            .map{ it -> it[0,1,2] }
                            .combine( fastq_paired_end_reads )
                            .combine( BOWTIE2_database.toList() )

    BOWTIE2_output = run_BOWTIE2_commands( BOWTIE2_run_params );

    //// BWA-MEM
    BWA_database = prepare_BWA_index ( params.amr_nucl_database ); 

    runs_ch.BWA_params
        .map{ it -> it[0,1,2] }
        .combine( fastq_paired_end_reads )
        .combine( BWA_database.toList() )
        .set{ BWA_run_params }

    BWA_output = run_BWA_commands( BWA_run_params );


    // GROOT
    // only get unique db parameters
    runs_ch.GROOT_params 
        .map{ it -> it[3] }
        .unique()
        .set{ GROOT_db_params_unique }

    GROOT_databases = prepare_GROOT_database( params.amr_nucl_database,
                                              GROOT_db_params_unique );
    runs_ch.GROOT_params 
        .combine( fastq_paired_end_reads )
        .combine( GROOT_databases.toList().toList() )
        .set{ GROOT_combined_run_params }

    GROOT_output = run_GROOT_commands( GROOT_combined_run_params );

    // ARIBA
    ARIBA_database = prepare_ARIBA_database( params.amr_database_version );

    runs_ch.ARIBA_params
        .map{ it -> it[0,1,2] }
        .combine( fastq_paired_end_reads )
        .combine( ARIBA_database )
        .set{ ARIBA_run_params }

    ARIBA_output = run_ARIBA_commands( ARIBA_run_params );

    // KMA
    runs_ch.KMA_params
        .map{ it -> it[1,3] }
        .set{ KMA_db_params }

    KMA_database = prepare_KMA_database( params.amr_nucl_database,
                                         KMA_db_params );

    runs_ch.KMA_params
        .combine( fastq_paired_end_reads )
        .combine( KMA_database.flatten().toList().toList() )
        .set{ KMA_combined_run_params }

    KMA_output = run_KMA_commands( KMA_combined_run_params );

    // HMMSEARCHNT
    HMMSEARCHNT_database = prepare_HMMSEARCHNT_database( params.amr_nucl_database ); 

    runs_ch.HMMSEARCHNT_params
        .map{ it -> it[0,1,2] }
        .combine( fasta_single_ended_reads )
        .combine( HMMSEARCHNT_database.toList() )
        .set{ HMMSEARCHNT_run_params }

    HMMSEARCHNT_output = run_HMMSEARCHNT_commands( HMMSEARCHNT_run_params );


    /*
    Nucleotide Methods vs Protein database:
            - BLASTX
            - DIAMONDBLASTX
            - PALADIN
            - MMSEQSBLASTX
    */


    // BLASTX
    BLASTX_database = prepare_BLASTX_BLASTP_database( params.amr_prot_database );

    // iterations of BLASTN params invocation to run correctly
    BLASTX_run_params = runs_ch.BLASTX_params
                                .map{ it -> it[0,1,2] }
                                .combine( fasta_single_ended_reads )
                                .combine( BLASTX_database.toList() )

    BLASTX_output = run_BLASTX_commands( BLASTX_run_params );
    
    // used for both DIAMONDBLASTX and DIAMONDBLASTP
    DIAMOND_database = prepare_DIAMOND_databases( params.amr_prot_database );

    runs_ch.DIAMONDBLASTX_params
        .map{ it -> it[0,1,2] }
        .combine( fastq_single_ended_reads )
        .combine( DIAMOND_database.toList() )
        .set{ DIAMONDBLASTX_run_params }

    DIAMONDBLASTX_output = run_DIAMONDBLASTX_commands( DIAMONDBLASTX_run_params );

    // PALADIN
    PALADIN_database = prepare_PALADIN_database( params.amr_prot_database );

    runs_ch.PALADIN_params
        .map{ it -> it[0,1,2] }
        .combine( fastq_single_ended_reads )
        .combine( PALADIN_database.toList() )
        .set{ PALADIN_run_params }

    PALADIN_output = run_PALADIN_commands( PALADIN_run_params, 
                                           params.amr_prot_database );

    // MMSEQSBLASTX
    MMSEQS_database = prepare_MMSEQS_database( params.amr_prot_database );

    runs_ch.MMSEQSBLASTX_params
        .map{ it -> it[0,1,2] }
        .combine( fastq_single_ended_reads )
        .combine( MMSEQS_database.db.toList() )
        .combine( MMSEQS_database.index.toList() )
        .set{ MMSEQSBLASTX_run_params }
    
    MMSEQSBLASTX_output = run_MMSEQSBLASTX_commands( MMSEQSBLASTX_run_params );
}
