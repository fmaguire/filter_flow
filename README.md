# AMR-metagenomics workflow

Workflow to search a metagenomic dataset for AMR genes using a large number of 
homology search tools and paramterisations. 

This is intended to work with the [CARD](https://card.mcmaster.ca/) database ([Alcock et al. 2020](https://www.ncbi.nlm.nih.gov/pubmed/31665441)) and only identifies homology related resistance (not mutational and/or ribosomal).

You must download a copy of the version you want to use from [here](https://card.mcmaster.ca/download)

It is possible to modify this to use another AMR database as long as there are nucleotide and protein sequence files provided.

## Usage

0. Requires a working install of [nextflow](https://www.nextflow.io/) and [conda](https://docs.conda.io/en/latest/).
All other dependencies are installed automatically.

1. To configure the workflow for your data, modify `params.config` accordingly:

- Provide paths to the paired fastq reads representing your metagenome as `reads_pe`

- Provide the read length (for groot) as `read_length`

- Provide the CARD version being used (for ARIBA) as `amr_database_version`

- Provide paths to the `nucleotide_fasta_protein_homolog_model.fasta` as `amr_nucl_database`

- Provide paths to the `protein_fasta_protein_homolog_model.fasta` as `amr_prot_database` 

e.g. 

```params {

    reads_pe = "$baseDir/metagenome/test_*{1,2}.fq"
    
    run_params_csv = "$baseDir/run_params.csv"
    
    amr_nucl_database = "$baseDir/database/nucleotide_fasta_protein_homolog_model.fasta"

    amr_prot_database = "$baseDir/database/protein_fasta_protein_homolog_model.fasta"
    
    amr_database_version = "3.0.9"

    read_length = 250

}```

2. You can modify `nextflow.config` to increase the number of CPUs being used and/or details of reported resource usage.

2. You can also modify `run_params.csv` as you see fit to add more/different paramterisations for tools being used.

3. Execute the workflow: `./nextflow metagenomic_amr_detection.nf -resume` (`-resume` will re-use cached results in the event of a crash and rerunning of the workflow).

## Outputs

Each tool output will be under the requisite folder in `results`:

```results/
├── aa
│   ├── blastp
│   ├── diamond-blastp
│   ├── hmm
│   └── mmseqsblastp
├── nt
│   ├── ariba
│   ├── blastn
│   ├── bowtie2
│   ├── groot
│   ├── hmm
│   └── kma
└── nt_to_aa
    ├── blastx
    ├── diamondblastx
    ├── mmseqsblastx
    └── paladin```
