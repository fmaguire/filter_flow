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
          prepare_HMMSEARCHNT_database; run_HMMSEARCHNT_commands;
        } from './modules/nt.nf'


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
          prepare_MMSEQS_database; run_MMSEQSBLASTX_commands;
        } from "./modules/nt_to_aa.nf"

/*
Protein Methods vs Protein database:
        - BLASTP
        - DIAMONDBLASTP
        - MMSEQSBLASTP
        - HMMSEARCHAA
*/

include { run_BLASTP_commands;  
          run_DIAMONDBLASTP_commands;
          run_MMSEQSBLASTP_commands; 
          prepare_HMMSEARCHAA_database; run_HMMSEARCHAA_commands; 
          } from "./modules/aa.nf"


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
                    DIAMONDBLASTX_params: it[0] == 'diamondblastx'
                    DIAMONDBLASTP_params: it[0] == 'diamondblastp'
                    KMA_params: it[0] == 'kma'
                    HMMSEARCHNT_params: it[0] == 'hmmsearchnt'
                    HMMSEARCHAA_params: it[0] == 'hmmsearchaa'
                    PALADIN_params: it[0] == 'paladin'
                    MMSEQSBLASTX_params: it[0] == 'mmseqsblastx'
                    MMSEQSBLASTP_params: it[0] == 'mmseqsblastp'
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

    BWA_run_params = runs_ch.BWA_params
                        .map{ it -> it[0,1,2] }
                        .combine( fastq_paired_end_reads )
                        .combine( BWA_database.toList() )

    BWA_output = run_BWA_commands( BWA_run_params );


    // GROOT
    // only get unique db parameters
    GROOT_db_params_unique = runs_ch.GROOT_params 
                                .map{ it -> it[3] }
                                .unique()

    GROOT_databases = prepare_GROOT_database( params.amr_nucl_database,
                                              GROOT_db_params_unique );
    GROOT_combined_run_params = runs_ch.GROOT_params 
                                        .combine( fastq_paired_end_reads )
                                        .combine( GROOT_databases.toList().toList() )

    GROOT_output = run_GROOT_commands( GROOT_combined_run_params );

    // ARIBA
    ARIBA_database = prepare_ARIBA_database( params.amr_database_version );

    ARIBA_run_params = runs_ch.ARIBA_params
                            .map{ it -> it[0,1,2] }
                            .combine( fastq_paired_end_reads )
                            .combine( ARIBA_database )

    ARIBA_output = run_ARIBA_commands( ARIBA_run_params );

    // KMA
    KMA_db_params = runs_ch.KMA_params
                        .map{ it -> it[1,3] }

    KMA_database = prepare_KMA_database( params.amr_nucl_database,
                                         KMA_db_params );

    KMA_combined_run_params = runs_ch.KMA_params
                                    .combine( fastq_paired_end_reads )
                                    .combine( KMA_database.flatten().toList().toList() )

    KMA_output = run_KMA_commands( KMA_combined_run_params );

    // HMMSEARCHNT
    HMMSEARCHNT_database = prepare_HMMSEARCHNT_database( params.amr_nucl_database ); 

    HMMSEARCHNT_run_params = runs_ch.HMMSEARCHNT_params
                                .map{ it -> it[0,1,2] }
                                .combine( fasta_single_ended_reads )
                                .combine( HMMSEARCHNT_database.toList() )

    HMMSEARCHNT_output = run_HMMSEARCHNT_commands( HMMSEARCHNT_run_params );


    /*
    Nucleotide Methods vs Protein database:
            - BLASTX
            - DIAMONDBLASTX
            - PALADIN
            - MMSEQSBLASTX
    */


    // BLASTX
    BLASTX_BLASTP_database = prepare_BLASTX_BLASTP_database( params.amr_prot_database );

    // iterations of BLASTN params invocation to run correctly
    BLASTX_run_params = runs_ch.BLASTX_params
                                .map{ it -> it[0,1,2] }
                                .combine( fasta_single_ended_reads )
                                .combine( BLASTX_BLASTP_database.toList() )

    BLASTX_output = run_BLASTX_commands( BLASTX_run_params );
    
    // used for both DIAMONDBLASTX and DIAMONDBLASTP
    DIAMOND_database = prepare_DIAMOND_databases( params.amr_prot_database );

    DIAMONDBLASTX_run_params = runs_ch.DIAMONDBLASTX_params
                                    .map{ it -> it[0,1,2] }
                                    .combine( fastq_single_ended_reads )
                                    .combine( DIAMOND_database.toList() )

    DIAMONDBLASTX_output = run_DIAMONDBLASTX_commands( DIAMONDBLASTX_run_params );

    // PALADIN
    PALADIN_database = prepare_PALADIN_database( params.amr_prot_database );

    PALADIN_run_params = runs_ch.PALADIN_params
                            .map{ it -> it[0,1,2] }
                            .combine( fastq_single_ended_reads )
                            .combine( PALADIN_database.toList() )

    PALADIN_output = run_PALADIN_commands( PALADIN_run_params, 
                                           params.amr_prot_database );

    // MMSEQSBLASTX
    MMSEQS_database = prepare_MMSEQS_database( params.amr_prot_database );


    MMSEQSBLASTX_run_params =runs_ch.MMSEQSBLASTX_params
        .map{ it -> it[0,1,2] }
        .combine( fastq_single_ended_reads )
        .combine( MMSEQS_database.db.toList() )
        .combine( MMSEQS_database.index.toList() )
    
    MMSEQSBLASTX_output = run_MMSEQSBLASTX_commands( MMSEQSBLASTX_run_params );

    /*
    Protein Methods vs Protein database:
        - BLASTP
        - DIAMONDBLASTP
        - MMSEQSBLASTP
        - HMMSEARCHAA
    */

    //BLASTP
    BLASTP_run_params = runs_ch.BLASTP_params
                            .map{ it -> it[0,1,2] }
                            .combine( fasta_single_ended_amino_acid )
                            .combine( BLASTX_BLASTP_database.toList() )

    BLASTP_output = run_BLASTP_commands( BLASTP_run_params );

    //DIAMONDBLASTP
    DIAMONDBLASTP_run_params = runs_ch.DIAMONDBLASTP_params
                                    .map{ it -> it[0,1,2] }
                                    .combine( fasta_single_ended_amino_acid )
                                    .combine( DIAMOND_database.toList() )

    DIAMONDBLASTP_output = run_DIAMONDBLASTP_commands( DIAMONDBLASTP_run_params );

    //MMSEQSBLASTP
    MMSEQSBLASTP_run_params = runs_ch.MMSEQSBLASTP_params
                                    .map{ it -> it[0,1,2] }
                                    .combine( fasta_single_ended_amino_acid )
                                    .combine( MMSEQS_database.db.toList() )
                                    .combine( MMSEQS_database.index.toList() )

    MMSEQSBLASTP_output = run_MMSEQSBLASTP_commands( MMSEQSBLASTP_run_params );

    //HMMSEARCHAA
    HMMSEARCHAA_database = prepare_HMMSEARCHAA_database( params.amr_prot_database ); 

    HMMSEARCHAA_run_params = runs_ch.HMMSEARCHAA_params
                                    .map{ it -> it[0,1,2] }
                                    .combine( fasta_single_ended_amino_acid )
                                    .combine( HMMSEARCHAA_database.toList() )

    HMMSEARCHAA_output = run_HMMSEARCHAA_commands( HMMSEARCHAA_run_params );

}
