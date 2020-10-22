/* 

    Parsers for each output into a unified format:
        read_name \t tool \t run_params \t aro_hit


*/

process parse_BLAST_output {
    
    input:
        path BLAST_output
    output:
        path "${BLAST_output}_parsed.tsv"
    
        # threshold the evalues 
    script:

        parser_BLAST_output

    path "${label}.out6" 

}

process parse_KMA_output {

    script:
        # add kmer size and threshold the mapping score
        """
        zcat output.frag.gz | awk -F $'\t' '{print $6, $7, $3}' 
        """

}
