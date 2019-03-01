#! /usr/bin/env python
'''
Identifies isoforms based on data in .wig files.
'''

import os

import numpy as np

import util


label = '_0.01'


def z_peak_score(wigs, half_width, gap, threshold):
    '''
    Calculates peak z scores for the wig data.
    Based on get_data_compute_statistic_20180130.m

    Args:
        wigs (ndarray[float]): 2D array, dims: (n strands, genome size), for reads on each strand
        half_width (int): half width of window for averaging for peak z score
        gap (int): gap left out (both sides) of central position for peak z score
        threshold (float): average read density (read/nt) threshold for consideration

    Returns:
        z_peak (ndarray[float]): 2D array, dims: (n strands, genome size) peak z score

    TODO:
        Add other filter options as in filter_generator_20180130.m?
    '''

    dims = wigs.shape
    z_peak_sum = np.zeros(dims)
    z_peak_prod = np.ones(dims)

    window = half_width * 2 + 1

    # Convolution vectors
    right_conv = np.zeros(window)
    left_conv = np.zeros(window)
    right_conv[half_width+1:] = 1 / half_width
    left_conv[:half_width] = 1 / half_width

    for i, strand in enumerate(wigs):
        strand2 = strand**2

        for conv in [right_conv, left_conv]:
            average = np.convolve(strand, conv, 'same')
            std = np.sqrt(np.convolve(strand2, conv, 'same') - average**2)
            value = (average > threshold) * (strand - average) / (std + (average == 0))
            z_peak_sum[i, :] += value
            value[value < 0] = 0
            z_peak_prod[i, :] *= value

    z_peak = 2 * z_peak_prod / (z_peak_sum + (z_peak_sum == 0))  # harmonic mean

    return z_peak

def z_step_score(wigs, half_width, gap, threshold):
    '''
    Calculates step z scores for the wig data.
    Based on get_data_compute_statistic_20180130.m

    Args:
        wigs (ndarray[float]): 2D array, dims: (n strands, genome size), for reads on each strand
        half_width (int): half width of window for averaging for step z score
        gap (int): gap left out (both sides) of central position for step z score
        threshold (float): average read density (read/nt) threshold for consideration

    Returns:
        z_step (ndarray[float]): 2D array, dims: (n strands, genome size) step z score

    TODO:
        Add other filter options as in filter_generator_20180130.m?
    '''

    dims = wigs.shape
    z_step = np.zeros(dims)

    window = (half_width + gap) * 2 + 1

    # Convolution vectors
    right_conv = np.zeros(window)
    left_conv = np.zeros(window)
    right_conv[half_width+2*gap+1:] = 1 / half_width
    left_conv[:half_width] = 1 / half_width

    for i, strand in enumerate(wigs):
        strand2 = strand**2

        right_average = np.convolve(strand, right_conv, 'same')
        left_average = np.convolve(strand, left_conv, 'same')
        right_std = np.sqrt(np.convolve(strand2, right_conv, 'same') - right_average**2)
        left_std = np.sqrt(np.convolve(strand2, left_conv, 'same') - left_average**2)
        positive_samples = (right_average > threshold) | (left_average > threshold)
        z_step[i, :] = positive_samples * (right_average - left_average) / np.sqrt(right_std + left_std + ~positive_samples)

    return z_step

def identify_regions(starts, ends, no_reads, gene_pad, ma_pad):
    '''
    Identify transcription unit region starts and ends based on no reads between regions.

    Args:
        starts (ndarray[int]): annotated gene start locations
        ends (ndarray[int]): annotated gene end locations
        no_reads (ndarray[int]): locations of no reads
        gene_pad (int): padding of base pairs to include before or after a gene
        ma_pad (int): pad around central point used for moving average

    Returns:
        real_starts (list[int]): location of starts of regions with reads (transcription units)
        real_ends (list[int]): location of ends of regions with reads (transcription units)
    '''

    real_starts = []
    real_ends = []

    n_genes = len(starts)
    n_no_reads = len(no_reads)
    gene = 0
    pos = 0

    while gene < n_genes:
        if gene == 0 and starts[gene] < no_reads[pos] - gene_pad:
            # Start at beginning if first gene starts before first period of no reads
            start = 1
        else:
            # Start at last point of no reads before gene otherwise
            start = no_reads[no_reads < starts[gene] - gene_pad][-1]
            pos = np.where(no_reads == start)[0][0]

        real_starts.append(max(0, start - gene_pad) + 1)

        while gene < n_genes - 1 and pos < n_no_reads - 1:
            if ends[gene] + gene_pad > no_reads[pos] - ma_pad:
                pos += 1
            elif no_reads[pos] + ma_pad > starts[gene + 1] - gene_pad:
                gene += 1
            else:
                real_ends.append(min(util.GENOME_SIZE, no_reads[pos] + gene_pad + 1))
                break

        gene += 1

    ## Add the end of the genome as the final end if not added already
    if len(real_starts) != len(real_ends):
        real_ends.append(util.GENOME_SIZE)

    return real_starts, real_ends


if __name__ == '__main__':
    # Parameters for analysis
    ## z statistics
    half_width_peak = 100
    half_width_step = 100
    threshold = 0.01
    gap_peak = 2  # not used unless filter is changed in z_score_statistics
    gap_step = 3

    ## Moving average of reads for identifying regions
    ma_pad = 5
    ma_window = 2*ma_pad + 1
    conv = np.ones(ma_window) / ma_window
    ma_min_reads = 0.5
    gene_pad = 100

    # Calculate Statistics
    wigs = util.load_wigs()
    total_reads = np.zeros((2, util.GENOME_SIZE))
    total_reads[0, :] = wigs[0, :] + wigs[2, :]
    total_reads[1, :] = wigs[1, :] + wigs[3, :]
    print('Calculating z scores...')
    z_peak = z_peak_score(wigs, half_width_peak, gap_peak, threshold)
    z_step = z_step_score(total_reads, half_width_step, gap_step, threshold)

    z_peak[z_peak < 1] = 1
    z_step_neg = -z_step
    z_step_pos = z_step
    z_step_neg[z_step_neg < 0.1] = 0.1
    z_step_pos[z_step_pos < 0.1] = 0.1

    # Load genome info
    genes, locus, all_starts, all_ends = util.load_genome()

    # Identify separate regions of genes for the forward strand
    print('Identifying forward regions...')
    reads = wigs[0, :] + wigs[2, :]
    ma_reads = np.convolve(reads, conv, 'same')
    no_reads = np.where(ma_reads < ma_min_reads)[0]

    starts = all_starts[all_starts > 0]
    ends = all_ends[all_ends > 0]

    real_starts, real_ends = identify_regions(starts, ends, no_reads, gene_pad, ma_pad)

    # Plot forward strand regions
    print('Plotting forward regions...')
    start_peak = z_peak[2, :]
    end_peak = z_peak[0, :]
    start_step = z_step_neg[0, :]
    end_step = z_step_pos[0, :]
    labels = ['z peak (start)', 'z peak (end)', 'z step (start)', 'z step (end)']
    for i, (start, end) in enumerate(zip(real_starts, real_ends)):
        scores = np.vstack((start_peak[start:end], end_peak[start:end], start_step[start:end], end_step[start:end]))
        util.plot_reads(start, end, genes, all_starts, all_ends, wigs, scores=scores, score_labels=labels,
            path=os.path.join(util.OUTPUT_DIR, f'fwd_{i}{label}.png'))

    # Identify separate regions of genes for the reverse strand
    print('Identifying reverse regions...')
    reads = wigs[1, :] + wigs[3, :]
    ma_reads = np.convolve(reads[::-1], conv, 'same')
    no_reads = np.where(ma_reads < ma_min_reads)[0]

    starts = -all_starts[all_starts < 0][::-1]
    ends = -all_ends[all_ends < 0][::-1]

    real_starts, real_ends = identify_regions(starts, ends, no_reads, gene_pad, ma_pad)

    # Plot reverse strand regions
    print('Plotting reverse regions...')
    start_peak = z_peak[3, ::-1]
    end_peak = z_peak[1, ::-1]
    start_step = z_step_pos[1, ::-1]
    end_step = z_step_neg[1, ::-1]
    labels = ['z peak (start)', 'z peak (end)', 'z step (start)', 'z step (end)']
    for i, (start, end) in enumerate(zip(real_starts[::-1], real_ends[::-1])):
        scores = np.vstack((start_peak[-end:-start], end_peak[-end:-start], start_step[-end:-start], end_step[-end:-start]))
        util.plot_reads(-start, -end, genes, all_starts, all_ends, wigs, scores=scores, score_labels=labels,
            path=os.path.join(util.OUTPUT_DIR, f'rev_{i}{label}.png'))
