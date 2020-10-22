/*

        Take input paired fastq reads and parse them into 
        per-tool channels.

        Secondarily, convert to single ended and fasta versions and load into
        channels for the tools that need input formatted that way.

*/

process build_single_ended {
    input:
        tuple val(label), path(reads)
    output:
        path "metagenome_se.fastq" 
    shell:
        """
        cat ${reads[0]} ${reads[1]} > metagenome_se.fastq
        """
}


process convert_to_fasta {
    input: 
        tuple val(label), path(reads)
    output:
        path "metagenome_reads.fasta" 
    shell:
        """
        cat ${reads[0]} ${reads[1]} | sed -n '1~4s/^@/>/p;2~4p' > metagenome_reads.fasta
        """
}


process convert_to_amino_acid {
    tag{ "ORFM: amino acid query generation" }
    conda "$baseDir/conda_envs/orfm.yml"
    input:
        path(reads_se) 
    output:
        path "metagenome_read_orfs.faa" 
    // 60 nucleotides => 30 amino acids
    shell:
        """
        orfm -m 60 ${reads_se} > metagenome_read_orfs.faa
        """
}

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


