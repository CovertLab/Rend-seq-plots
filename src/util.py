'''
Useful functions and constants for scripts.
'''

import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
WIGS_DIR = os.path.join(DATA_DIR, 'wigs')
GENOME_DIR = os.path.join(DATA_DIR, 'genome')
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'out')
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

# .wig files from Lalanne et al
WIG_FILE_TEMPLATE = os.path.join(WIGS_DIR, 'GSM2971252_Escherichia_coli_WT_Rend_seq_5_exo_MOPS_comp_25s_frag_pooled_{}_no_shadow.wig')
PROCESSED_WIG_FILE = os.path.join(DATA_DIR, 'GSM2971252_reads.npy')
WIG_STRANDS = ['3f', '3r', '5f', '5r']

# Genome information
GENOME_SIZE = 4639675
ANNOTATION_FILE = os.path.join(GENOME_DIR, 'U00096.2.faa')


def load_wigs(cached=None, wig_file=None, strands=None, genome_size=None):
    '''
    Loads read data from .wig files

    Args:
        cached (str): path to cached data, if None, defaults to PROCESSED_WIG_FILE.
            If file exists, it will be read, otherwise data will be loaded from
            the raw data and written to file.
        wig_file (str): template for the wig files to be read (gets formatted for each strand),
            if None, defaults to WIG_FILE_TEMPLATE
        strands (list[str]): each strand name that gets inserted into the wig_file string
            (order determines order of first dimension of returned array), if None, defaults
            to WIG_STRANDS
        genome_size (int): base pair length of the genome, if None, defaults to GENOME_SIZE

    Returns:
        ndarray[float]: 2D array, dims: (n strands, genome size), for reads on each strand
    '''

    # Check for cached results and return if exists
    if cached is None:
        cached = PROCESSED_WIG_FILE
    if os.path.exists(cached):
        print('Loaded cached read data')
        return np.load(cached)

    # Defaults in file for E. coli if none given
    if wig_file is None:
        wig_file = WIG_FILE_TEMPLATE
    if strands is None:
        strands = WIG_STRANDS
    if genome_size is None:
        genome_size = GENOME_SIZE

    data = np.zeros((len(strands), genome_size))
    for i, strand in enumerate(strands):
        with open(wig_file.format(strand), 'r') as f:
            reader = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONNUMERIC)

            # skip lines of text at top of the file
            d = []
            while len(d) == 0:
                try:
                    d = np.array(list(reader))
                except Exception as exc:
                    print(exc)

            index = np.array(d[:, 0], dtype=int) - 1
            reads = np.array(d[:, 1], dtype=float)

            data[i, index] = reads

    np.save(cached, data)

    return data

def load_genome(path=None):
    '''
    Loads information about the genome, specifically genes and their locations.

    Args:
        path (str): path to genome file, if None, defaults to ANNOTATION_FILE

    Returns:
        genes (ndarray[str]): gene names
        locus_tags (ndarray[str]): locus tag names
        starts (ndarray[int]): genome position of start of gene
        ends (ndarray[int]): genome position of end of gene
    '''

    # Set defaults
    if path is None:
        path = ANNOTATION_FILE

    # Read in file
    with open(path) as f:
        data = f.read()

    # Initialize data structures to return
    genes = []
    locus_tags = []
    starts = []
    ends = []

    # Read in each gene and pull out relevant information
    for gene in data.split('\n>'):
        locus = re.findall('\[locus_tag=(.*?)\]', gene)[0]
        gene_name = re.findall('\[gene=(.*?)\]', gene)[0]

        location = re.findall('\[location=(.*?)\]', gene)[0]
        # Handle complement strand
        if location[0] == 'c':
            dir = -1
        else:
            dir = 1
        start, end = (dir * int(pos) for pos in re.findall('([0-9]*)\.\.([0-9]*)', location)[0])

        genes.append(gene_name)
        locus_tags.append(locus)
        starts.append(start)
        ends.append(end)

    genes = np.array(genes)
    locus_tags = np.array(locus_tags)
    starts = np.array(starts)
    ends = np.array(ends)

    sort_idx = np.argsort(starts)

    return genes[sort_idx], locus_tags[sort_idx], starts[sort_idx], ends[sort_idx]

def z_score_statistics(wigs, half_width_peak=100, half_width_step=100, threshold=0.25, gap_peak=2, gap_step=3):
    '''
    Calculates peak and step z scores for the wig data.
    Based on get_data_compute_statistic_20180130.m

    Args:
        wigs (ndarray[float]): 2D array, dims: (n strands, genome size), for reads on each strand
        half_width_peak (int): half width of window for averaging for peak z score
        half_width_step (int): half width of window for averaging for step z score
        threshold (float): average read density (read/nt) threshold for consideration
        gap_peak (int): gap left out (both sides) of central position for peak z score
        gap_step (int): gap left out (both sides) of central position for step z score

    Returns:
        z_peak (ndarray[float]): 2D array, dims: (n strands, genome size) peak z score
        z_step (ndarray[float]): 2D array, dims: (n strands, genome size) step z score

    TODO:
        Add other filter options as in filter_generator_20180130.m?
    '''

    dims = wigs.shape

    z_peak = np.ones(dims)
    z_step = np.zeros(dims)

    window_peak = half_width_peak * 2 + 1
    window_step = (half_width_step + gap_step) * 2 + 1

    # Convolution vectors for peak calculation
    right_conv_peak = np.zeros(window_peak)
    left_conv_peak = np.zeros(window_peak)
    right_conv_peak[half_width_peak+1:] = 1 / half_width_peak
    left_conv_peak[:half_width_peak] = 1 / half_width_peak

    # Convolution vectors for step calculation
    right_conv_step = np.zeros(window_step)
    left_conv_step = np.zeros(window_step)
    right_conv_step[half_width_step+2*gap_step+1:] = 1 / half_width_step
    left_conv_step[:half_width_step] = 1 / half_width_step

    for i, strand in enumerate(wigs):
        strand2 = strand**2

        # Calculate peak z score
        for conv in [right_conv_peak, left_conv_peak]:
            average = np.convolve(strand, conv, 'same')
            std = np.sqrt(np.convolve(strand2, conv, 'same') - average**2)
            z_peak[i, :] *= (average > threshold) * (strand - average) / (std + (average == 0))

        # Calculate step z score
        right_average = np.convolve(strand, right_conv_step, 'same')
        left_average = np.convolve(strand, left_conv_step, 'same')
        right_std = np.sqrt(np.convolve(strand2, right_conv_step, 'same') - right_average**2)
        left_std = np.sqrt(np.convolve(strand2, left_conv_step, 'same') - left_average**2)
        positive_samples = (right_average > threshold) | (left_average > threshold)
        z_step[i, :] = np.abs(positive_samples * (right_average - left_average) / np.sqrt(right_std + left_std + ~positive_samples))

    return z_peak, z_step


def plot_reads(start, end, genes, starts, ends, reads, fit=None, path=None):
    '''
    Plots the reads of the 3' and 5' data on the given strand.  Also shows any
    genes that start or finish within the specified region.

    Args:
        start (int): start position in the genome
        end (int): end position in the genome
        genes (array of str): names of genes
        starts (array of int): start position for each gene
        ends (array of int): end position for each gene
        reads (2D array of float): reads for each strand at each position
            dims (strands x genome length)
        path (str): path to save image, if None, just displays image to screen
    '''

    # Set forward (0) or reverse (1) strand to correspond to data order
    if start < 0:
        strand = 1
    else:
        strand = 0

    # Identify genes in region and get reads
    # Need to subtract 1 from index for reads since starts at 0 not 1 like genome position
    if strand:
        # Reverse strand needs indices adjusted
        mask = (-ends > -start) & (-starts < -end)
        loc = np.arange(end, start)
        with np.errstate(divide='ignore'):
            three_prime = np.log(reads[strand, int(-start-1):int(-end-1)][::-1])
            five_prime = np.log(reads[2+strand, int(-start-1):int(-end-1)][::-1])
    else:
        # Forward strand
        mask = (ends > start) & (starts < end)
        loc = np.arange(start, end)
        with np.errstate(divide='ignore'):
            three_prime = np.log(reads[strand, int(start-1):int(end-1)])
            five_prime = np.log(reads[2+strand, int(start-1):int(end-1)])

    genes = genes[mask]
    starts = starts[mask]
    ends = ends[mask]

    # Adjust for 0 reads from log
    three_prime[three_prime < 0] = -1
    five_prime[five_prime < 0] = -1

    plt.figure()
    plt.step(loc, np.vstack((three_prime, five_prime)).T, linewidth=0.25)

    if fit is not None:
        with np.errstate(divide='ignore'):
            plt.step(loc, np.log(fit), color='k')

    gene_line = -1.5
    gene_offset = 0.1
    for gene, s, e in zip(genes, starts, ends):
        plt.plot([s, e], [gene_line, gene_line], 'k')
        plt.plot([s, s], [gene_line-gene_offset, gene_line+gene_offset], 'k')
        plt.plot([e, e], [gene_line-gene_offset, gene_line+gene_offset], 'k')
        plt.text((s+e)/2, gene_line-3*gene_offset, gene, ha='center', fontsize=6)
    plt.xlim([loc[0], loc[-1]])

    plt.xlabel('Genome Location')
    plt.ylabel('Reads (log)')

    if path is None:
        plt.show()
    else:
        plt.savefig(path)
    plt.close('all')
