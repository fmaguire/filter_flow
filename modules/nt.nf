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
        path(amr_ref)
    output:
        path 'amr.blastn.db.*'
    script:
        """
        makeblastdb -dbtype nucl -in $amr_ref -out amr.blastn.db
        """
}


//BLASTN defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/
process run_BLASTN_commands {
    tag { "BLASTN: ${label}" }
    publishDir "results/nt/blastn", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/blast.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(read_fasta), path(blastn_database) 
    output:
        path "${label}.out6" 
    script:
        """
        blastn -max_target_seqs 1 -query $read_fasta $run_params -num_threads ${task.cpus} -db amr.blastn.db -outfmt 6 > ${label}.out6
        """
}

// BOWTIE2
process prepare_BOWTIE2_index {
    conda "$baseDir/conda_envs/bowtie2.yml"
    input:
        path(amr_ref)
    output:
        path 'amr.bowtie2.db*' 
    script:
        """
        bowtie2-build --threads ${task.cpus} ${amr_ref} amr.bowtie2.db
        """
}

process run_BOWTIE2_commands {
    tag { "BOWTIE2: ${label}" }
    publishDir "results/nt/bowtie2", pattern: "${label}.sam", mode: 'copy'
    conda "$baseDir/conda_envs/bowtie2.yml"
    input:
        tuple val(tool), val(label), val(run_params), val(read_label), path(reads), path(bowtie2_index)
    output:
        path("${label}.sam")
    script:
        """
        bowtie2 -x amr.bowtie2.db $run_params -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} > ${label}.sam
        """
}


// BWA-MEM
process prepare_BWA_index {
    conda "$baseDir/conda_envs/bwa.yml"
    input:
        path amr_ref 
    output:
        path 'amr.bwa.db*' 
    script:
        """
        bwa index -p amr.bwa.db ${amr_ref} 
        """
}


process run_BWA_commands {
    tag { "BWA-MEM: ${label}" }
    publishDir "results/nt/bwa", pattern: "${label}.sam", mode: "copy"
    conda "$baseDir/conda_envs/bwa.yml"
    input:
        tuple val(tool), val(label), val(run_params), val(read_label), path(reads), path(bwa_index) 
    output:
        path "${label}.sam" 
    script:
        """
        bwa mem -t ${task.cpus} $run_params amr.bwa.db ${reads[0]} ${reads[1]} > ${label}.sam
        """
}
// GROOT

process prepare_GROOT_database {
    tag { "!GROOT_db: ${db_params}" }
    conda "$baseDir/conda_envs/groot.yml"
    input:
        path amr_ref 
        val db_params 
    output:
        path 'groot_db_*' 
    script:
        """
        mkdir clustered_amr_db;
        vsearch --cluster_size $amr_ref --id 0.90 --msaout MSA.tmp;
        awk '!a[\$0]++ {of="clustered_amr_db/cluster-" ++fc ".msa"; print \$0 >> of ; close(of)}' RS= ORS="\\n\\n" MSA.tmp && rm MSA.tmp;
        groot index ${db_params} -m clustered_amr_db -i groot_db_${db_params[-2,-1]} -w ${params.read_length};
        """
}

process run_GROOT_commands {
    conda "$baseDir/conda_envs/groot.yml"
    tag { "GROOT: ${label}" }
    publishDir "results/nt/groot", pattern: "${label}.bam", mode: "copy"

    input:
        tuple val(tool), val(label), val(run_params), val(db_params), val(read_label), path(reads), path(groot_db)
    output:
        path "${label}.bam" 
    script:
        """
        groot align ${run_params} -i groot_db_${db_params[-2,-1]} -p ${task.cpus} -f ${reads[0]} ${reads[1]} > ${label}.bam
        """
}

// ARIBA 
process prepare_ARIBA_database {
    conda "$baseDir/conda_envs/ariba.yml"
    input:
        val amr_db_version 
    output:
        path 'ariba_amr_db' 
    // test at awk live generation for increase portability to other databases in future
    // cat ${amr_ref} | awk -F'|' '/^>/ {split(\$6,name," "); split(\$5,aro,":"); output_str= sprintf("%s.%s.%s.%s.%s\\t1\\t0\\t.\\t.\\t%s", name[1], aro[2], $2, $4, NR, name[1]); print output_str}' > ariba_meta.tsv
    script:
        """
        ariba getref --version ${amr_db_version} card card_database
		ariba prepareref -f card_database.fa -m card_database.tsv ariba_amr_db
        """
}

process run_ARIBA_commands {
    conda "$baseDir/conda_envs/ariba.yml"
    tag{ "ARIBA: ${label}" }
	publishDir "results/nt/ariba", pattern: "${label}_output", mode: "copy"
	input:
        tuple val(tool), val(label), val(run_params), val(read_label), path(reads), path(ariba_db)
	output:	
		path "${label}_output" 
	script:
		"""
		ariba run --noclean ${run_params} ${ariba_db} ${reads[0]} ${reads[1]} ${label}_output 
		"""
}

// KMA
process prepare_KMA_database {
    conda "$baseDir/conda_envs/kma.yml"
    input:
        path amr_ref 
        tuple val(label), val(db_params)
    output:
        path "kma_database_${label}.*"
    script:
        """
        kma index -i ${amr_ref} ${db_params} -o kma_database_${label}
        """
}

process run_KMA_commands {
    conda "$baseDir/conda_envs/kma.yml"
    tag{ "KMA: ${label}" }
    publishDir "results/nt/kma", pattern: "${label}.frag.gz", mode: "copy"
    input:
        tuple val(tool), val(label), val(run_params), val(db_params), val(read_label), path(reads), path(kma_db)
    output:
        path "${label}.frag.gz"
    script:
        """
        kma -ipe ${reads[0]} ${reads[1]} -t ${task.cpus} -1t1 -o ${label} -t_db kma_database_${label}
        """
}

// HMMSearch nucleotide
process prepare_HMMSEARCHNT_database {
    conda "$baseDir/conda_envs/hmm.yml"
	publishDir "results/hmm_cluster_nt", pattern: "*.msa_clean", mode: "copy"
    input:
        path amr_ref 
    output:
        path 'card90.hmm', emit: hmm
		path ("*.msa_clean"), emit: clusters
    script:
        """
        vsearch --cluster_size ${amr_ref} --id 0.90 --msaout MSA.tmp
        awk '!a[\$0]++ {of="./cluster-" ++fc ".msa"; print \$0 >> of ; close(of)}' RS= ORS="\\n\\n" MSA.tmp && rm MSA.tmp;
        for msa in *.msa;
            do 
                sed 's/+/-/g' \${msa} > \${msa}_clean
                hmmbuild --fragthresh 1 \${msa}.hmm \${msa}_clean
            done
        cat *.hmm > card90.hmm
        """
}

process run_HMMSEARCHNT_commands {
    conda "$baseDir/conda_envs/hmm.yml"
    tag{ "HMM: ${label}" }
    publishDir "results/nt/hmm", pattern: "${label}.tbl", mode: "copy"
    input:
        tuple val(tool), val(label), val(run_params), path(read_fasta), path(hmm_db) 
    output:
        path "${label}.tbl"
    script: 
        """
        hmmsearch --cpu ${task.cpus} --noali --notextw --tblout ${label}.tbl ${hmm_db} ${read_fasta}
        """
}

