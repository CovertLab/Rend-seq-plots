#! /usr/bin/env python
"""
Plots R-end seq read counts for a specified region in the E. coli genome.
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

import src.util as util


def plot_read_counts(range, reversed, wig_index, ymax, filename):
    """
    Plots R-end seq read counts for a specified region in the E. coli genome.
    """
    reads = util.load_wigs(wig_index=wig_index)

    start = range[0]
    end = range[1]

    # Reverse strand
    if reversed:
        loc = np.arange(end, start)
        three_prime = reads[1, (-start - 1):(-end - 1)][::-1]
        five_prime = reads[3, (-start - 1):(-end - 1)][::-1]
    # Forward strand
    else:
        loc = np.arange(start, end)
        three_prime = reads[0, (start - 1):(end - 1)]
        five_prime = reads[2, (start - 1):(end - 1)]

    genes, _, gene_starts, gene_ends = util.load_genome()
    mask = (gene_starts >= start) & (gene_ends <= end)

    fig, axes = plt.subplots(2, 1, figsize=(9, 3), sharex=True, gridspec_kw={'height_ratios': [7, 1]})

    # Plot reads
    ax0 = axes[0]
    ax0.step(loc, three_prime, linewidth=1, alpha=0.5, zorder=1, label="3'-mapped")
    ax0.step(loc, five_prime, linewidth=1, alpha=0.5, zorder=0, label="5'-mapped")
    ax0.set_xlim([start, end])
    ax0.set_ylim([1, 10**ymax])
    ax0.set_ylabel('R-end seq reads')
    ax0.set_yscale('log')
    ax0.axes.get_xaxis().set_visible(False)
    ax0.spines['left'].set_position(('outward', 15))
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.legend(loc=1)

    # Plot gene locations
    ax1 = axes[1]
    for gene, gene_start, gene_end in zip(genes[mask], gene_starts[mask], gene_ends[mask]):
        ax1.hlines(0, gene_start, gene_end, linewidth=2)
        ax1.vlines([gene_start, gene_end], -1, 1, linewidth=1)
        ax1.text((gene_start + gene_end)/2, 2, gene, size=8, ha='center')

    ax1.set_ylim([-3, 3])
    ax1.axes.get_yaxis().set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.ticklabel_format(axis='x', style='plain')
    ax1.set_xlabel('Genome position (nt)')

    plt.tight_layout()
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    plt.savefig(os.path.join(out_dir, filename))


def parse_args():
    """
    Parses arguments from the command line.
    Returns:
        ArgumentParser namespace: values of variables parsed from the command line
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--range',
        nargs=2,
        type=int,
        metavar=('start', 'end'),
        help='Range of genomic coordinates to plot')
    parser.add_argument('-w', '--wig',
        type=int,
        default=0,
        help=f'Index of wig type to use (0-{len(util.WIG_FILE_TEMPLATES)-1}, default: 0)')
    parser.add_argument('--ymax',
        type=int,
        default=5,
        help='Max value of y axis in log10 scale')
    parser.add_argument('-f', '--filename',
        type=str,
        default='output.png',
        help='Name of output file')
    parser.add_argument('--reversed',
        dest='rev',
        action='store_true',
        help='If set, plot reads in the reversed direction')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    plot_read_counts(args.range, args.rev, args.wig, args.ymax, args.filename)
