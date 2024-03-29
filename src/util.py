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
ECOCYC_DIR = os.path.join(DATA_DIR, 'ecocyc')
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'out')
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

# .wig files from Lalanne et al
WIG_FILE_TEMPLATES = [
    os.path.join(WIGS_DIR, 'GSM2500131_Escherichia_coli_WT_Rend_seq_MOPS_comp_25s_frag_pooled_{}_no_shadow.wig'),
    os.path.join(WIGS_DIR, 'GSM2971252_Escherichia_coli_WT_Rend_seq_5_exo_MOPS_comp_25s_frag_pooled_{}_no_shadow.wig'),
    os.path.join(WIGS_DIR, 'GSM2971253_Escherichia_coli_pnp_Rend_seq_MOPS_comp_25s_frag_{}.wig'),
    os.path.join(WIGS_DIR, 'GSM2971254_Escherichia_coli_rnb_Rend_seq_MOPS_comp_25s_frag_{}.wig'),
    ]
PROCESSED_WIG_FILES = [
    os.path.join(DATA_DIR, 'GSM2500131_reads.npy'),
    os.path.join(DATA_DIR, 'GSM2971252_reads.npy'),
    os.path.join(DATA_DIR, 'GSM2971253_reads.npy'),
    os.path.join(DATA_DIR, 'GSM2971254_reads.npy'),
    ]
WIG_STRANDS = ['3f', '3r', '5f', '5r']

# Genome information
GENOME_SIZE = 4639675
ANNOTATION_FILE = os.path.join(GENOME_DIR, 'U00096.2.faa')


def load_wigs(wig_index=1, strands=None, genome_size=None, cached=True):
    '''
    Loads read data from .wig files

    Args:
        wig_index (int): index for wig template to use from WIG_FILE_TEMPLATES
        strands (list[str]): each strand name that gets inserted into the wig_file string
            (order determines order of first dimension of returned array), if None, defaults
            to WIG_STRANDS
        genome_size (int): base pair length of the genome, if None, defaults to GENOME_SIZE
        cached (bool): if True, loads cached data if it exists

    Returns:
        ndarray[float]: 2D array, dims: (n strands, genome size), for reads on each strand
    '''

    # Check for cached results and return if exists
    cached_file = PROCESSED_WIG_FILES[wig_index]
    if cached and os.path.exists(cached_file):
        return np.load(cached_file)

    # Defaults in file for E. coli if none given
    if strands is None:
        strands = WIG_STRANDS
    if genome_size is None:
        genome_size = GENOME_SIZE

    wig_file = WIG_FILE_TEMPLATES[wig_index]
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

    np.save(cached_file, data)

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

def load_annotated_tu_directions():
    '''
    Gets annotated transcription unit directions (fwd or rev) from ecocyc data.
    Fwd strand is represented by 1, rev strand is represented by -1.

    Returns:
        dict {TU (str): direction (int)}: annotated direction for a given transcription unit
    '''

    tu_file = os.path.join(ECOCYC_DIR, 'tus.tsv')

    tu_directions = {}
    with open(tu_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip headers
        for line in reader:
            if '&larr' in line[1]:
                direction = -1
            elif '&rarr' in line[1]:
                direction = 1
            else:
                continue

            tu_directions[line[0]] = direction

    return tu_directions

def load_annotated_promoter_positions():
    '''
    Gets annotated promoter start positions from ecocyc data.

    Returns:
        ndarray[int]: annotated promoter positions
    '''

    promoter_file = os.path.join(ECOCYC_DIR, 'promoters.tsv')
    tus = load_annotated_tu_directions()

    starts = []
    with open(promoter_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip headers
        for line in reader:
            pos = line[1]
            info = line[2]
            if pos:
                for tu in info.split(' // '):
                    if tu in tus:
                        starts.append(int(pos) * tus[tu])

    return np.unique(starts)

def load_annotated_terminator_positions():
    '''
    Gets annotated terminator positions from ecocyc data.

    Returns:
        ndarray[int]: annotated terminator positions
    '''

    terminator_file = os.path.join(ECOCYC_DIR, 'terminators.tsv')
    tus = load_annotated_tu_directions()

    terms = []
    with open(terminator_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # skip headers
        for line in reader:
            left_pos = line[1]
            right_pos = line[2]
            info = line[3]
            for tu in info.split(' // '):
                if tu in tus:
                    direction = tus[tu]
                    if direction == 1 and right_pos:
                        terms.append(int(right_pos))
                    elif direction == -1 and left_pos:
                        terms.append(-int(left_pos))

    return np.unique(terms)

def fit_annotations_to_data(annotations, reads, window=10, read_cutoff=4):
    '''
    Adjust annotated positions to better match spikes in seq data by finding the max read counts
    in a window around the annotated position.

    Args:
        annotations (ndarray[int]): annotated positions for promoter or terminators
        reads (ndarray[float]): 2D array for fwd and rev reads, dims (2 (fwd, rev) x GENOME_SIZE)
            5' wig data (5f or 5r) for promoters
            3' wig data (3f or 3r) for terminators
        winodw (int): window to either side of point to check

    Returns:
        ndarray[int]: adjusted annotations to align with data peaks
            positive values will be forward strand starting with 1
            negative values will be reverse strand starting with -1
    '''

    adjusted_annotations = []

    for pos in annotations:
        # Reverse strand
        if pos < 0:
            start = -(pos + window + 1)
            end = -(pos - window)
            data = reads[1, start:end]

            max_pos = np.argmax(data)
            if data[max_pos] < read_cutoff:
                continue

            adjusted_annotations.append(-(start + max_pos + 1))
        # Forward strand
        else:
            start = pos - window - 1
            end = pos + window
            data = reads[0, start:end]

            max_pos = np.argmax(data)
            if data[max_pos] < read_cutoff:
                continue

            adjusted_annotations.append(start + max_pos + 1)

    return np.unique(adjusted_annotations)

def plot_reads(start, end, genes, starts, ends, reads, scores=None, score_labels=None, threshold=5, path=None):
    '''
    Plots the reads of the 3' and 5' data on the given strand.  Also shows any
    genes that start or finish within the specified region.

    Args:
        start (int): start position in the genome
        end (int): end position in the genome
        genes (ndarray[str]): names of genes
        starts (ndarray[int]): start position for each gene
        ends (ndarray[int]): end position for each gene
        reads (ndarray[float]): reads for each strand at each position
            dims (strands x genome length)
        scores (ndarray[float]): statistic scores for each position in the region
            to be plotted, dims (n scores x region length), if None, only reads are plotted
        score_labels (list[str]): label for each score given, if None, no labels shown
        threshold (float): threshold level to indicate on plots
        path (str): path to save image, if None, just displays image to screen
    '''

    def plot_genes(plt, genes, starts, ends):
        gene_line = -1.5
        gene_offset = 0.1
        for gene, s, e in zip(genes, starts, ends):
            plt.plot([s, e], [gene_line, gene_line], 'k')
            plt.plot([s, s], [gene_line - gene_offset, gene_line + gene_offset], 'k')
            plt.plot([e, e], [gene_line - gene_offset, gene_line + gene_offset], 'k')
            plt.text((s + e) / 2, gene_line - 3 * gene_offset, gene, ha='center', fontsize=6)
        plt.xlim([loc[0], loc[-1]])

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

    if scores is not None:
        n_scores = scores.shape[0]
    else:
        n_scores = 0
    n_plots = n_scores + 1

    plt.figure(figsize=(8.5, 4*n_plots))

    # Plot only reads
    ax = plt.subplot(n_plots, 1, 1)
    ax.step(loc, np.vstack((three_prime, five_prime)).T, linewidth=0.25)
    plot_genes(plt, genes, starts, ends)
    plt.ylabel('Reads (log)')

    # Plot each score on a separate subplot
    if scores is not None:
        for i, score in enumerate(scores):
            score[score < 0.1] = 0.1
            ax = plt.subplot(n_plots, 1, i+2)
            ax.step(loc, np.vstack((three_prime, five_prime)).T, linewidth=0.25)
            with np.errstate(divide='ignore'):
                ax.step(loc, np.log(score), color='k')

            if score_labels is not None:
                plt.ylabel(score_labels[i])
            plot_genes(plt, genes, starts, ends)
            plt.axhline(np.log(threshold), color='r', linestyle='--')

    plt.xlabel('Genome Location')

    if path is None:
        plt.show()
    else:
        plt.savefig(path)
    plt.close('all')

def plot_tus(start, end, genes, starts, ends, reads, tus, annotated, path=None):
    '''
    Plots the reads of the 3' and 5' data on the given strand and the transcription
    units identified.

    Args:
        start (int): start position in the genome
        end (int): end position in the genome
        genes (ndarray[str]): names of genes
        starts (ndarray[int]): start position for each gene
        ends (ndarray[int]): end position for each gene
        reads (ndarray[float]): reads for each strand at each position
            dims (strands x genome length)
        tus (list[set[str]]): set of all unique transcription units in each gene region
            with genes in the same transcription unit separated by :
        path (str): path to save image, if None, just displays image to screen
    '''

    def plot_genes(plt, genes, starts, ends):
        gene_line = -1.5
        gene_offset = 0.1
        for gene, s, e in zip(genes, starts, ends):
            plt.plot([s, e], [gene_line, gene_line], 'k')
            plt.plot([s, s], [gene_line - gene_offset, gene_line + gene_offset], 'k')
            plt.plot([e, e], [gene_line - gene_offset, gene_line + gene_offset], 'k')
            plt.text((s + e) / 2, gene_line - 3 * gene_offset, gene, ha='center', fontsize=6)
        plt.xlim([loc[0], loc[-1]])

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
    plot_genes(plt, genes, starts, ends)

    if annotated:
        p = load_annotated_promoter_positions()
        t = load_annotated_terminator_positions()

        if strand:
            for pos in p[(p < start) & (p > end)]:
                plt.axvline(pos, color='g', linestyle=(0, (4, 6)))

            for pos in t[(t < start) & (t > end)]:
                plt.axvline(pos, color='r', linestyle=(0, (4, 6)))
        else:
            for pos in p[(p > start) & (p < end)]:
                plt.axvline(pos, color='g', linestyle=(0, (4, 6)))

            for pos in t[(t > start) & (t < end)]:
                plt.axvline(pos, color='r', linestyle=(0, (4, 6)))


    tu_start = -2
    tu_step = 0.2
    for i, tu in enumerate(sorted(tus)):
        idx = np.array([np.where(genes == g)[0][0] for g in tu.split(':')])
        start = np.min(starts[idx])
        end = np.max(ends[idx])

        pos = tu_start - tu_step * i
        plt.plot([start, end], [pos, pos], 'r', linewidth=2)

    plt.ylabel('Reads (log)')
    plt.xlabel('Genome Location')

    if path is None:
        plt.show()
    else:
        plt.savefig(path)
    plt.close('all')
