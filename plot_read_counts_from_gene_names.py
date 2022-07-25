#! /usr/bin/env python
"""
Plots R-end seq read counts from a .tsv file containing gene names.
"""
import argparse
import numpy as np
import os

import src.util as util
from plot_read_counts import plot_read_counts


def plot_read_counts_from_gene_names(input_file, wig_index, ymax):
    """
    Plots R-end seq read counts from a .tsv file containing gene names.
    """
    genes, _, gene_starts, gene_ends = util.load_genome()
    gene_name_to_start = {name: start for (name, start) in zip(genes, gene_starts)}
    gene_name_to_end = {name: end for (name, end) in zip(genes, gene_ends)}

    with open(input_file, 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            operon_genes = line[:-1].split('\t')[:-1]
            try:
                operon_gene_starts = [gene_name_to_start[gene] for gene in operon_genes]
                operon_gene_ends = [gene_name_to_end[gene] for gene in operon_genes]
            except KeyError:
                print(f'Could not find Rend-seq data for operon {"-".join(operon_genes)}.')
                for gene in operon_genes:
                    if gene not in gene_name_to_start:
                        print(gene)
                continue

            rev = operon_gene_starts[0] > operon_gene_starts[1]
            if rev:
                range = [operon_gene_starts[0] + 100, operon_gene_ends[1] - 100]
                filename = '-'.join(operon_genes[::-1]) + '_w3.pdf'
            else:
                range = [operon_gene_starts[0] - 100, operon_gene_ends[1] + 100]
                filename = '-'.join(operon_genes) + '_w3.pdf'

            plot_read_counts(range, rev, wig_index, ymax, filename)



def parse_args():
    """
    Parses arguments from the command line.
    Returns:
        ArgumentParser namespace: values of variables parsed from the command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
        type=str,
        default='data/wcEcoli/operon_table.tsv',
        help='Name of input file containing gene names')
    parser.add_argument('-w', '--wig',
        type=int,
        default=0,
        help=f'Index of wig type to use (0-{len(util.WIG_FILE_TEMPLATES)-1}, default: 0)')
    parser.add_argument('--ymax',
        type=int,
        default=5,
        help='Max value of y axis in log10 scale')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    plot_read_counts_from_gene_names(args.input, args.wig, args.ymax)
