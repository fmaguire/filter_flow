/*
Nucleotide Methods vs Protein database:
        - BLASTX
        - DIAMONDBLASTX
        - PALADIN
        - MMSEQSBLASTX
*/


// BLASTX
process prepare_BLASTX_BLASTP_database {
    conda "$baseDir/conda_envs/blast.yml"
    input:
        path amr_ref 
    output:
        path 'amr.prot.db.*' 
    script:
        """
        makeblastdb -dbtype prot -in $amr_ref -out amr.prot.db
        """
}

//BLASTX defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options/
process run_BLASTX_commands {
    tag { "BLASTX: ${label}" }
    publishDir "results/nt_to_aa/blastx", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/blast.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(read_fasta), path(blastx_database)
    output:
        path "${label}.out6" 
    script:
        """
        blastx -max_target_seqs 1 -query $read_fasta $run_params -num_threads ${task.cpus} -db amr.prot.db -outfmt 6 > ${label}.out6
        """
} 

// DIAMOND_BLASTX
process prepare_DIAMOND_databases {
    conda "$baseDir/conda_envs/diamond.yml"
    input:
        path amr_ref 
    output:
        path 'amr.diamond.db.*' 
    script:
        """
        diamond makedb --in $amr_ref --db amr.diamond.db
        """
}

//DIAMOND BLASTX defaults: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options/
process run_DIAMONDBLASTX_commands {
    tag { "DIAMONDBLASTX: ${label}" }
    publishDir "results/nt_to_aa/diamondblastx", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/diamond.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(reads_se), path(diamond_database)
    output:
        path "${label}.out6" 
    script:
        """
        diamond blastx -p ${task.cpus} --max-target-seqs 1 --db ${diamond_database} --outfmt 6 --out ${label}.out6 ${run_params} -q ${reads_se}
        """
}

// PALADIN
process prepare_PALADIN_database {
    conda "$baseDir/conda_envs/paladin.yml"
    input:
        path amr_ref
    output:
        path "${amr_ref}.*"
    script:
        """
        paladin index -r3 ${amr_ref} 
        """
}

process run_PALADIN_commands {
    tag { "PALDIN: ${label}" }
    publishDir "results/nt_to_aa/paladin", pattern: "${label}.bam", mode: "copy"
    conda "$baseDir/conda_envs/paladin.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(reads_se), path(paladin_database)
        path(amr_ref)
    output:
        path "${label}.bam" 
    script:
        """
        paladin align -F 0.5 -t ${task.cpus} -a -C -V ${run_params} ${amr_ref} ${reads_se} | samtools view -Sb > ${label}.bam
        """
}

// MMSEQSBLASTX

process prepare_MMSEQS_database {
    conda "$baseDir/conda_envs/mmseqs.yml"
    input:
        path amr_ref 
    output:
        path 'mmseqs_db*', emit: db
        path 'mmseqs_tmp', emit: index
    script:
        """
        mmseqs createdb ${amr_ref} mmseqs_db 
        mmseqs createindex mmseqs_db mmseqs_tmp
        """
}

process run_MMSEQSBLASTX_commands {
    tag { "MMSEQSBLASTX: ${label}" }
    publishDir "results/nt_to_aa/mmseqsblastx", pattern: "${label}.out6", mode: "copy"
    conda "$baseDir/conda_envs/mmseqs.yml"
    input:
        tuple val(tool), val(label), val(run_params), path(reads_se), path(mmseqs_database), path(mmseqs_index)
    output: 
        path "${label}.out6"
    script:
        """
        mmseqs easy-search --format-mode 2 --threads ${task.cpus} ${run_params} ${reads_se} mmseqs_db ${label}.out6 ${mmseqs_index}
        """
}

