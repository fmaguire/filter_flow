#!/usr/bin/env python

import subprocess
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os
import numpy as np
from scipy.special import logsumexp
import shutil
from pathlib import Path
import collections

def parse_args():
    """
    Parser user CLI args
    """
    parser = argparse.ArgumentParser(description="Simulate a metagenome from"
                                                 "genomes, RGI annotations, "
                                                 "and AMR genes to insert")
    parser.add_argument("-m", "--metadata", required=True,
                         help="TSV containing 3 columns: path to genome, "
                              "path to RGI annotation for genome, and desired"
                              "relative abundance for genome (no header)")

    parser.add_argument("-a", "--amr_insert",
                         help="TSV containing AMR sequences to insert (last 3 "
                              "columns optional, without them insertion "
                              "will be random): 'Name Accession Sequence "
                              "(Contig Location Strand)' (no header) "
                              " if you provide location check it isn't within "
                              " AMR gene already in contig")

    parser.add_argument("-r", "--read_error_profile",
                        help="ART error profile to use",
                        default="MSv3",
                        choices=['MSv3', 'HS25'])

    parser.add_argument('--frag_size',
                        type=int,
                        help="Fragment size for paired end reads",
                        default=500)

    parser.add_argument('--frag_std',
                        type=float,
                        help="Fragment size standard deviation for paired end "
                             "reads",
                        default=50)

    parser.add_argument('--random_seed',
                        type=int,
                        help="Random seed to use for simulation and insertion",
                        default=42)

    parser.add_argument("-o", "--output_dir",
                        default="output",
                        help="Directory to output results to")

    parser.add_argument("-l", "--read_length",
                        type=int,
                        default=250,
                        help="Read length to simulate")

    parser.add_argument("--minimum_overlap",
                        type=int,
                        default=50,
                        help="Minimum read overlap to label as AMR gene")

    parser.add_argument("--num_reads",
                        type=int,
                        default=25000000,
                        help="Number of reads to simulate")

    args = parser.parse_args()

    return args


def prepare_data(metadata_fp):
    """
    Parse provided metadata TSV, parse the linked rgi output and copy
    genomes for AMR insertions
    """
    metadata = pd.read_csv(metadata_fp, sep='\t', names=['genome_fp',
                                                         'rgi_fp',
                                                         'relative_abundance'])

    # check if all the filepaths can be found
    genome_paths_exist = metadata['genome_fp'].apply(lambda x: \
                                                        os.path.exists(x))
    if not genome_paths_exist.all():
        missing = metadata.loc[~genome_paths_exist, 'genome_fp'].values
        raise ValueError(f"Genomes can't be found: {missing}")

    rgi_annotation_exist = metadata['rgi_fp'].apply(lambda x: \
                                                        os.path.exists(x))
    if not genome_paths_exist.all():
        missing = metadata.loc[~rgi_annotation_exist, 'rgi_fp'].values
        raise ValueError(f"RGI annotations can't be found: {missing}")

    # dictionary of existing annotations across genomes keyed by contig
    # and containing start, end, amr_name, and ARO for each seq
    amr_annotations = collections.defaultdict(list)
    for rgi_output in metadata['rgi_fp']:
        rgi_data = pd.read_csv(rgi_output, sep='\t')
        for _, row in rgi_data.iterrows():
            contig = row['Contig']
            data = {'start': row['Start'], 'stop': row['Stop'],
                    'aro': row['ARO'], 'amr_name': row['Best_Hit_ARO'],
                    'strand': row['Orientation']}
            amr_annotations[contig].append(data)

    # contig names must be unique
    genome_to_contig = collections.defaultdict(list)
    sequence_data = {}
    for genome_fp in metadata['genome_fp']:
        genome = Path(genome_fp).name
        for record in SeqIO.parse(genome_fp, 'fasta'):
            genome_to_contig[genome].append(record.id)

            if record.id in sequence_data:
                raise ValueError(f"{record.id} is duplicated, all contig name"
                                  " must be unique")
            else:
                sequence_data[record.id] = record.seq

    return metadata, genome_to_contig, sequence_data, amr_annotations


def insert_amr_genes(amr_insert_fp, sequence_data, amr_annotations):
    """
    Check number of fields in AMR_insert_data and randomly generate them if
    they are not defined
    """
    amr_insert_data = pd.read_csv(amr_insert_fp, sep='\t',
                                  names=['name', 'accession',
                                         'nt_sequence', 'contig',
                                         'location', 'strand'])
    all_contigs = list(sequence_data.keys())
    contig_lengths = {contig: len(sequence_data[contig]) for \
                                                contig in sequence_data}

    # randomly select contig, location, and strand if not defined
    # could be implemented using apply but iteration is clearer and should
    # not be too slow
    for ix, row in amr_insert_data.iterrows():
        if pd.isna(row['contig']):
            contig = np.random.choice(all_contigs)
            amr_insert_data.loc[ix, 'contig'] = contig

        # forbid insertion into existing AMR gene
        if pd.isna(row['location']):
            contig_length = contig_lengths[contig]
            existing_amr_mask = np.zeros((contig_length,))
            for annotation in amr_annotations[contig]:
                if annotation['start'] < annotation['stop']:
                    existing_amr_mask[slice(annotation['start'],
                                            annotations['stop'])] = -100
                else:
                    existing_amr_mask[slice(annotation['stop'],
                                            annotations['start'])] = -100

            # gumbel-max might be faster but this is actually intelligible
            # in a year or two, blame Dr. G. Gray (pers. comm.)
            # select uniformly from non-masked/blastlisted sites
            prob = np.exp(existing_amr_mask - logsumexp(existing_amr_mask))
            insert_loc = np.random.choice(np.arange(contig_length),
                                          p=prob)
            amr_insert_data.loc[ix, 'location'] = int(insert_loc)

        if pd.isna(row['strand']):
            amr_insert_data.loc[ix, 'strand'] = np.random.choice(['+', '-'])

    # insert and adjust the amr_annotations accordingly
    for ix, row in amr_insert_data.iterrows():
        contig = row['contig']
        insert_loc = int(row['location'])
        insert_strand = row['strand']
        if insert_strand == '-':
            insert_seq = Seq(row['nt_sequence']).reverse_complement()
            insert_length = len(insert_seq)
            insert_start = insert_loc + insert_length
            insert_stop = insert_loc
        elif insert_strand == '+':
            insert_seq = Seq(row['nt_sequence'])
            insert_length = len(insert_seq)
            insert_start = insert_loc
            insert_stop = insert_loc + insert_length
        else:
            raise ValueError(f"Invalid strand for {row} in amr_insert_data")

        # find all existing genes in contig
        for gene in amr_annotations[contig]:
            # if existing gene starts after AMR gene then add length to index
            if gene['start'] > insert_loc and gene['stop'] > insert_loc:
                gene['start'] = gene['start'] + insert_length
                gene['stop'] = gene['stop'] + insert_length
            elif gene['start'] < insert_loc and gene['stop'] < insert_loc:
                # as the insert comes after we don't have to do anything
                pass
            else:
                raise ValueError(f"AMR gene being inserted in {contig} within"
                                     f" {gene} at loc: {row['location']}")

        # add amr gene to sequence_data
        contig_seq = sequence_data[contig]
        contig_seq = contig_seq[:insert_loc] + insert_seq + contig_seq[insert_loc:]
        sequence_data[contig] = contig_seq

        insert_data = {'start': insert_start,
                       'stop': insert_stop,
                       'strand': insert_strand,
                       'aro': row['accession'],
                       'amr_name': row['name']}
        amr_annotations[contig].append(insert_data)

    return sequence_data, amr_annotations


def duplicate_to_relative_abundance(metadata, genome_to_contig, sequence_data,
                                    amr_annotations, output_dir):
    """
    Based on relative abundance copy contigs the appropriate amount
    """

    # for relative abundance in metadata find all the contigs:
    for ix, row in metadata.iterrows():
        contigs = genome_to_contig[

    #


if __name__ == '__main__':

    args = parse_args()

    # parse all the already existing AMR annotations
    metadata, genome_to_contig, \
            sequence_data, amr_annotations = prepare_data(args.metadata)

    # if amr genes need inserted then add them
    if args.amr_insert:
        sequence_data, amr_annotations = insert_amr_genes(args.amr_insert,
                                                          sequence_data,
                                                          amr_annotations)
    # duplicate contigs as per relative abundance
    duplicate_to_relative_abundance(metadata, genome_to_contig, sequence_data,
                                    amr_annotations)



