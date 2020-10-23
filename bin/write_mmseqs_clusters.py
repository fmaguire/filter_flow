#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO

def run():
    parser = argparse.ArgumentParser(description='Write mmseq sequence cluster files and dump singeltons.')
    parser.add_argument('-c', '--cluster_fasta', type=str, required=True,
                    help="Path to MMSEQS2 easy-cluster fasta")
    parser.add_argument('-o', '--output_folder', type=str, required=True,
                    help="Folder to output family sequences")
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.mkdir(args.output_folder)

    # go inefficient to avoid bugs
    all_clusters = []
    cluster = []
    for record in SeqIO.parse(args.cluster_fasta, 'fasta'):
        if len(record.seq) == 0 and len(cluster) > 0:
            all_clusters.append(cluster)
            cluster = []
        else:
            cluster.append(record)

    clusterid = 0

    for cluster in all_clusters:
        output_fp = os.path.join(args.output_folder,
                                 str(clusterid) + ".faa")

        with open(output_fp, 'w') as out_fh:
            if len(cluster) == 1:
                cluster = cluster * 2
                cluster[0].id = cluster[0].id + "_copy"
                SeqIO.write(cluster, out_fh, 'fasta')
            else:
                SeqIO.write(cluster, out_fh, 'fasta')
        clusterid += 1

if __name__ == '__main__':
    run()

