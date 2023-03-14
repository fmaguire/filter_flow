/*
Protein Methods vs Protein database:
        - BLASTP
        - DIAMONDBLASTP
        - MMSEQSBLASTP
        - HMMSEARCHAA
*/

// * BLASTP
//BLASTP defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastp_application_options/
process run_BLASTP_commands {
    tag { "BLASTP: ${label}" }
    publishDir "results/aa/blastp", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/blast.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(read_aa_fasta), path(blast_database) 
    output:
        path "${label}.out6"
    script:
        """
        blastp -max_target_seqs 1  -query $read_aa_fasta $run_params -num_threads ${task.cpus} -db amr.prot.db -outfmt 6 > ${label}.out6
        """
}

// DIAMONDBLASTP
process run_DIAMONDBLASTP_commands {
    tag { "DIAMONDBLASTP: ${label}" }
    publishDir "results/aa/diamond-blastp", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/diamond.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(read_aa_fasta), path(diamond_database) 
    output:
        path "${label}.out6" 
    script:
        """
        diamond blastp -p ${task.cpus} --max-target-seqs 1 --db ${diamond_database} --outfmt 6 --out ${label}.out6 ${run_params} -q ${read_aa_fasta}
        """
} 

//MMSEQSBLASTP
process run_MMSEQSBLASTP_commands {
    tag { "MMSEQSBLASTP: ${label}" }
    publishDir "results/aa/mmseqsblastp", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/mmseqs.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(reads_aa_fasta), path(mmseqs_database), path(mmseqs_index)
    output: 
        path "${label}.out6"
    script:
        """
        mmseqs easy-search --format-mode 2 --threads ${task.cpus} ${run_params} ${reads_aa_fasta} mmseqs_db ${label}.out6 ${mmseqs_index}
        """
}

// HMMSEARCHAA
process prepare_HMMSEARCHAA_database {
    conda "$baseDir/conda_envs/hmm.yml"
	publishDir "results/hmm_cluster_aa", pattern: "cluster_dir/*.faa", mode: "copy"
    input:
        path amr_ref 
    output:
        path 'card90.hmm', emit: hmm
		path ("cluster_dir/*"), emit: clusters
    script:
        """
        mmseqs easy-cluster --min-seq-id 0.7 --remove-tmp-files 1 --threads ${task.cpus} ${amr_ref} clusters tmp
        mkdir -p cluster_dir
        write_mmseqs_clusters.py -c clusters_all_seqs.fasta -o cluster_dir
        for seq in cluster_dir/*.faa;
            do 
                mafft --auto \$seq > \${seq}.msa
                hmmbuild --fragthresh 1.0 \${seq}.hmm \${seq}.msa 
        done
        cat cluster_dir/*.hmm > card90.hmm
        """
}

process run_HMMSEARCHAA_commands {
    conda "$baseDir/conda_envs/hmm.yml"
    tag{ "HMM: ${label}" }
    publishDir "results/aa/hmm", pattern: "${label}.tbl", mode: "copy"
    input:
        tuple val(tool), val(label), val(run_params), path(read_aa_fasta), path(hmm_db) 
    output:
        path "${label}.tbl"
    script: 
        """
        hmmsearch --cpu ${task.cpus} --noali --notextw --tblout ${label}.tbl ${hmm_db} ${read_aa_fasta}
        """
}

