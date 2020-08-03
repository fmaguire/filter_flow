#!/usr/bin/env nextflow

/*

        Take input paired fastq reads and parse them into 
        per-tool channels.

        Secondarily, convert to single ended and fasta versions and load into
        channels for the tools that need input formatted that way.

*/

//// download CARD database
//process get_CARD {
//	input:
//        val card_canonical_version from params.card_version
//	output:
//		path "card/nucleotide_fasta_protein_homolog_model.fasta" into amr_nucl
//		path "card/protein_fasta_protein_homolog_model.fasta" into amr_prot
//	script:
//		"""
//        wget -P canonical_${card_canonical_version} https://card.mcmaster.ca/download/0/broadstreet-v{params.canonical_version}.tar.bz2 2>&1 >> {log}
//        tar -C canonical_${card_canonical_version} -xvf canonical_${card_canonical_version}/broadstreet-v${card_canonical_version}.tar.bz2",
//		"""
//}


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
metagenome_fastq_se.into{ aminoacid_conversion;
                         DIAMONDBLASTX_input}

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


// convert input metagenome to amino acid orfs using ORFM
process convert_to_amino_acid {
    tag{ "ORFM: amino acid query generation" }
    conda "$baseDir/conda_envs/orfm.yml"
    input:
        path(reads_se) from aminoacid_conversion
    output:
        path "metagenome_read_orfs.faa" into metagenome_aminoacids
    // 60 nucleotides => 30 amino acids
    shell:
        """
        orfm -m 60 ${reads_se} > metagenome_read_orfs.faa
        """
}
metagenome_aminoacids.into { BLASTP_input;
                             DIAMONDBLASTP_input;
                             HMMSEARCHPROTEIN_input;
                             MMSEQSPROTEIN_input }

/* 
        Parse run parameter files for sweep settings to use in 
        homology search 

*/
Channel
    .fromPath(params.run_params_csv)
    .splitCsv(header:true)
    .map{ row -> tuple(row.tool, row.label, row.run_params, row.db_params) }
    .branch {
        BLASTN_params: it[0] == "blastn"
        BLASTX_params: it[0] == "blastx"
        BOWTIE2_params: it[0] == "bowtie2"
        BWA_params: it[0] == "bwa"
        GROOT_params: it[0] == "groot"
		ARIBA_params: it[0] == "ariba"
        DIAMONDBLASTX_params: it[0] == 'diamond-blastx'
        DIAMONDBLASTP_params: it[0] == 'diamond-blastp'
     }
    .set{ runs_ch }


/*
Nucleotide Methods vs Nucleotide database:
        - BLASTN
        - bowtie2
        - BWA-MEM
        - groot
        - ariba

*/

// BLASTN
process prepare_BLASTN_database {
    conda "$baseDir/conda_envs/blast.yml"
    input:
        path amr_ref from params.amr_nucl_database
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
    .map{ it -> it[0,1,2] }
    .combine( BLASTN_input )
    .combine( BLASTN_database.toList() )
    .set{ BLASTN_run_params }

//BLASTN defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/
process run_BLASTN_commands {
    tag { "BLASTN: ${label}" }
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
        path amr_ref from params.amr_nucl_database
    output:
        path 'amr.bowtie2.db*' into BOWTIE2_database
    script:
        """
        bowtie2-build --threads ${task.cpus} ${amr_ref} amr.bowtie2.db
        """
}

// filter to just bowtie2
runs_ch.BOWTIE2_params
    .map{ it -> it[0,1,2] }
    .combine( BOWTIE2_input )
    .combine( BOWTIE2_database.toList() )
    .set{ BOWTIE2_run_params }

process run_bowtie2_commands {
    tag { "BOWTIE2: ${label}" }
    publishDir "results/nt/bowtie2", pattern: '*.bam', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/bowtie2.yml"
    input:
        set val(tool), val(label), val(run_params), val(read_label), path(reads), path(bowtie2_index) from BOWTIE2_run_params
    output:
        file '*.sam' into BOWTIE2_output
    script:
        """
        bowtie2 -x amr.bowtie2.db $run_params -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} > ${label}.sam
        """
}


// BWA-MEM
process prepare_bwa_index {
    conda "$baseDir/conda_envs/bwa.yml"
    input:
        path amr_ref from params.amr_nucl_database
    output:
        path 'amr.bwa.db*' into BWA_database
    script:
        """
        bwa index -p amr.bwa.db ${amr_ref} 
        """
}

runs_ch.BWA_params
    .map{ it -> it[0,1,2] }
    .combine( BWA_input )
    .combine( BWA_database.toList() )
    .set{ BWA_run_params }

process run_bwa_commands {
    tag { "BWA-MEM: ${label}" }
    publishDir "results/nt/bowtie2", pattern: '*.sam', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/bwa.yml"
    input:
        set val(tool), val(label), val(run_params), val(read_label), path(reads), path(bwa_index) from BWA_run_params
    output:
        file '*.sam' into BWA_output
    script:
        """
        bwa mem -t ${task.cpus} $run_params amr.bwa.db ${reads[0]} ${reads[1]} > ${label}.sam
        """
}

// GROOT
runs_ch.GROOT_params
    .into{ GROOT_db_params; GROOT_run_params }

// only get unique db parameters
GROOT_db_params
    .map{ it -> it[3] }
    .unique()
    .set{ GROOT_db_params_unique }

process prepare_groot_database {
    tag { "!GROOT_db: ${db_params}" }
    conda "$baseDir/conda_envs/groot.yml"
    input:
        path amr_ref from params.amr_nucl_database
        val db_params from GROOT_db_params_unique
    output:
        path 'groot_db_*' into GROOT_databases
    script:
        """
        mkdir clustered_amr_db;
        vsearch --cluster_size $amr_ref --id 0.90 --msaout MSA.tmp;
        awk '!a[\$0]++ {of="clustered_amr_db/cluster-" ++fc ".msa"; print \$0 >> of ; close(of)}' RS= ORS="\\n\\n" MSA.tmp && rm MSA.tmp;
        groot index ${db_params} -m clustered_amr_db -i groot_db_${db_params[-2,-1]} -w ${params.read_length};
        """
}

GROOT_run_params
    .combine( GROOT_input )
    .combine( GROOT_databases.toList().toList() )
    .set{ GROOT_combined_run_params }

process run_groot_commands {
    conda "$baseDir/conda_envs/groot.yml"
    tag { "GROOT: ${label}" }
    publishDir "results/nt/groot", pattern: '*.bam', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}

    input:
        set val(tool), val(label), val(run_params), val(db_params), val(read_label), path(reads), path(groot_db) from GROOT_combined_run_params
    output:
        file "*.bam" into GROOT_output
    script:
        """
        groot align ${run_params} -i groot_db_${db_params[-2,-1]} -p ${task.cpus} -f ${reads[0]} ${reads[1]} > ${label}.bam
        """
}

// ARIBA 
process prepare_ariba_database {
    //conda "$baseDir/conda_envs/ariba.yml"
    input:
        val amr_db_version from params.amr_database_version
    output:
        path 'ariba_amr_db' into ARIBA_database

    // test at awk live generation for increase portability to other databases in future
    // cat ${amr_ref} | awk -F'|' '/^>/ {split(\$6,name," "); split(\$5,aro,":"); output_str= sprintf("%s.%s.%s.%s.%s\\t1\\t0\\t.\\t.\\t%s", name[1], aro[2], $2, $4, NR, name[1]); print output_str}' > ariba_meta.tsv

    script:
        """
        ariba getref --version ${amr_db_version} card card_database
		ariba prepareref -f card_database.fa -m card_database.tsv ariba_amr_db
        """
}

runs_ch.ARIBA_params
    .map{ it -> it[0,1,2] }
    .combine( ARIBA_input )
    .combine( ARIBA_database )
    .set{ ARIBA_run_params }


process run_ariba_comamnds {
    //conda "$baseDir/conda_envs/ariba.yml"
    tag{ "ARIBA: ${label}" }
	publishDir "results/nt/ariba", pattern: "*_output"
	input:
        set val(tool), val(label), val(run_params), val(read_label), path(reads), path(ariba_db) from ARIBA_run_params
	output:	
		path "*_output" into ARIBA_output
	script:
		"""
		ariba run --noclean ${run_params} ${ariba_db} ${reads[0]} ${reads[1]} ${label}_output 
		"""
}

//// HMMSearch nucleotide
//process prepare_hmmsearch_nt_database {
//    conda "$baseDir/conda_envs/hmmsearch.yml"
//    input:
//        path amr_ref from params.amr_nucl_database
//        val db_params from GROOT_db_params_unique
//    output:
//        path 'groot_db_*' into GROOT_databases
// 
//
//}
//
//


/*
Nucleotide Methods vs Protein database:
        - BLASTX
        - DIAMONDBLASTX
        - PALADIN
*/


// BLASTX
process prepare_BLASTX_database {
    conda "$baseDir/conda_envs/blast.yml"
    input:
        path amr_ref from params.amr_prot_database
    output:
        path 'amr.blastx.db.*' into BLASTX_database
    script:
        """
        makeblastdb -dbtype prot -in $amr_ref -out amr.blastx.db
        """
}

// filter params to just blastx and combine params with input/db to get all
// iterations of BLASTX params invocation to run correctly
runs_ch.BLASTX_params
    .map{ it -> it[0,1,2] }
    .combine( BLASTX_input )
    .combine( BLASTX_database.toList() )
    .set{ BLASTX_run_params }

//BLASTX defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options/
process run_BLASTX_commands {
    tag { "BLASTX: ${label}" }
    publishDir "results/nt_to_aa/blastx", pattern: '*.out6', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/blast.yml"
    input:
        set val(tool), val(label), val(run_params), path(read_fasta), path(blastx_database) from BLASTX_run_params
    output:
        file '*.out6' into BLASTX_output
    script:
        """
        blastx -query $read_fasta $run_params -num_threads ${task.cpus} -db amr.blastx.db -outfmt 6 > ${label}.out6
        """
} 


// DIAMONDBLASTX
process prepare_DIAMOND_databases {
    conda "$baseDir/conda_envs/diamond.yml"
    input:
        path amr_ref from params.amr_prot_database
    output:
        path 'amr.diamond.db.*' into DIAMONDBLASTX_database, DIAMONDBLASTP_database
    script:
        """
        diamond makedb --in $amr_ref --db amr.diamond.db
        """
}

// filter params to just blastx and combine params with input/db to get all
// iterations of BLASTX params invocation to run correctly
runs_ch.DIAMONDBLASTX_params
    .map{ it -> it[0,1,2] }
    .combine( DIAMONDBLASTX_input )
    .combine( DIAMONDBLASTX_database.toList() )
    .set{ DIAMONDBLASTX_run_params }

//DIAMOND BLASTX defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options/
process run_DIAMONDBLASTX_commands {
    tag { "DIAMONDBLASTX: ${label}" }
    publishDir "results/nt_to_aa/diamond-blastx", pattern: '*.out6', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/diamond.yml"
    input:
        set val(tool), val(label), val(run_params), path(reads_se), path(diamond_database) from DIAMONDBLASTX_run_params
    output:
        file '*.out6' into DIAMONDBLASTX_output
    script:
        """
        diamond blastx -p ${task.cpus} --max-target-seqs 1 --db ${diamond_database} --outfmt 6 --out ${label}.out6 ${run_params} -q ${reads_se}
        """
} 


/*
Protein Methods vs Protein database:
        - BLASTP
        - DIAMONDBLASTP
        - MMSEQS
*/


runs_ch.DIAMONDBLASTP_params
    .map{ it -> it[0,1,2] }
    .combine( DIAMONDBLASTP_input )
    .combine( DIAMONDBLASTP_database.toList() )
    .set{ DIAMONDBLASTP_run_params }

//DIAMOND BLASTX defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options/
process run_DIAMONDBLASTP_commands {
    tag { "DIAMONDBLASTP: ${label}" }
    publishDir "results/aa/diamond-blastp", pattern: '*.out6', saveAs: { "${file(it).getSimpleName()}.${file(it).getExtension()}"}
    conda "$baseDir/conda_envs/diamond.yml"
    input:
        set val(tool), val(label), val(run_params), path(read_aa_fasta), path(diamond_database) from DIAMONDBLASTP_run_params
    output:
        file '*.out6' into DIAMONDBLASTP_output
    script:
        """
        diamond blastp -p ${task.cpus} --max-target-seqs 1 --db ${diamond_database} --outfmt 6 --out ${label}.out6 ${run_params} -q ${read_aa_fasta}
        """
} 
