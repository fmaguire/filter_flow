/*

        Takes the output from each tool and create summaries of ARGs detected with truth appended.
*/

process summarize_results {
    conda "$baseDir/conda_envs/postprocessing.yml"
	publishDir "results/", pattern: "ResultSummary/*.tsv", mode: "copy"
	input:
		path BLASTN_output
		path BOWTIE2_output
		path BWA_output
		path GROOT_output
		path KMA_output
		path HMMSEARCHNT_output
		path BLASTX_output
		path DIAMONDBLASTX_output
		path PALADIN_output
		path MMSEQSBLASTX_output
		path BLASTP_output
		path DIAMONDBLASTP_output
		path MMSEQSBLASTP_output
		path HMMSEARCHAA_output
		//^these input dont actually matter, it's just here to halt this process until everything else finishes
		
		path card_protein
		path card_dna
		path card_index
		path metagenome_labels
	output:
		path ("ResultSummary/*.tsv"), emit: results
    shell:
        """
		shopt -s extglob
		find $baseDir/results/!(hmm_cluster*) -type f -follow -print > files.txt
		shopt -u extglob
        assign_labels.py -f "files.txt" -d ${card_dna} -p ${card_protein} -i ${card_index} -t ${metagenome_labels} -hmmdna "$baseDir/results/hmm_cluster_nt" -hmmprot "$baseDir/results/hmm_cluster_aa/cluster_dir" -o ResultSummary 

        """
}