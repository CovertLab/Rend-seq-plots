#! /usr/bin/env python
'''
Identifies isoforms based on data in .wig files.
'''

import argparse
import os

import numpy as np

import util


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
        region_starts (list[int]): location of starts of regions with reads (transcription units)
        region_ends (list[int]): location of ends of regions with reads (transcription units)
    '''

    region_starts = []
    region_ends = []

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

        region_starts.append(max(0, start - gene_pad) + 1)

        while gene < n_genes - 1 and pos < n_no_reads - 1:
            if ends[gene] + gene_pad > no_reads[pos] - ma_pad:
                pos += 1
            elif no_reads[pos] + ma_pad > starts[gene + 1] - gene_pad:
                gene += 1
            else:
                region_ends.append(min(util.GENOME_SIZE, no_reads[pos] + gene_pad + 1))
                break

        gene += 1

    ## Add the end of the genome as the final end if not added already
    if len(region_starts) != len(region_ends):
        region_ends.append(util.GENOME_SIZE)

    return region_starts, region_ends

def identify_tu_locations(starts, ends, start_score, end_score, threshold=5, rev=False):
    '''
    Identifies transcription units from scores.

    Args:
        starts (list[int]): start position for each region
        ends (list[int]): end position for each region
        start_score (ndarray[float]): score for transcript starts at each position
        end_score (ndarray[float]): score for transcript ends at each position
        threshold (float): threshold for score to call a peak
        rev (bool): True if for reverse strand

    Returns:
        tu_starts (list[list[int]]): start position of each transcription unit for
            each gene region
        tu_ends (list[list[int]]): end position of each transcription unit for
            each gene region
    '''

    def combine_positions(score, threshold, gap=3):
        '''
        Combines nearby positions that are above the threshold to prevent duplicate
        positions for the same peak (start/stop).

        Args:
            score (ndarray[float]): score for start or end within a region of interest
            threshold (float): threshold for the score to be labeled a peak
            gap (int): maximum gap of positions being below threshold to consider all
                positions to be part of the same peak

        Returns:
            ndarray[int]: unique peak positions
        '''

        combined = []
        pos = np.where(score > threshold)[0]

        if len(pos) > 0:
            group = [pos[0]]
            for p in pos[1:]:
                if p - group[-1] < gap:
                    group.append(p)
                else:
                    combined.append(group[np.argmax(score[group])])
                    group = [p]

            combined.append(group[np.argmax(score[group])])
        return np.array(combined)

    all_starts = []
    all_ends = []

    for s, e in zip(starts, ends):
        threshold_starts = s - 1 + combine_positions(start_score[s-1:e-1], threshold)
        threshold_ends = s - 1 + combine_positions(end_score[s-1:e-1], threshold)

        tu_starts = []
        tu_ends = []
        for ts in threshold_starts:
            if rev:
                ts = -ts
            for te in threshold_ends:
                if rev:
                    te = -te

                if ts < te:
                    tu_starts.append(ts)
                    tu_ends.append(te)

        all_starts.append(tu_starts)
        all_ends.append(tu_ends)

    return all_starts, all_ends

def identify_tu_genes(tu_starts, tu_ends, genes, gene_starts, gene_ends):
    '''
    Identifies genes that are in a transcription unit together.

    Args:
        tu_starts (list[list[int]]): start position of each transcription unit for
            each gene region
        tu_ends (list[list[int]]): end position of each transcription unit for
            each gene region
        genes (ndarray[str]): gene names
        gene_starts (ndarray[int]): genome position of start of gene
        gene_ends (ndarray[int]): genome position of end of gene

    Returns:
        list[set[str]]: set of all unique transcription units in each gene region
            with genes in the same transcription unit separated by :
    '''

    tu_genes = []

    for starts, ends in zip(tu_starts, tu_ends):
        tus = set()
        for start, end in zip(starts, ends):
            tu = genes[(gene_starts >= start) & (gene_starts <= end) & (gene_ends >= start) & (gene_ends <= end)]
            if len(tu) > 0:
                tus.add(':'.join(tu))

        tu_genes.append(tus)

    return tu_genes

def parse_args():
    '''
	Parses arguments from the command line.

	Returns:
		ArgumentParser namespace: values of variables parsed from the command line
    '''

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--only-fwd',
        dest='rev',
        action='store_false',
        help='Only run analysis on the forward strand')
    parser.add_argument('-r', '--only-rev',
        dest='fwd',
        action='store_false',
        help='Only run analysis on the reverse strand')
    parser.add_argument('--no-plot',
        dest='plot',
        action='store_false',
        help='Skip plotting figures')
    parser.add_argument('-z', '--plot-z',
        action='store_true',
        help='Plot z statistic plots')
    parser.add_argument('-l', '--label',
        default='',
        type=str,
        help='Label to append to plot outputs')
    parser.add_argument('-g', '--gene',
        help='Search for region and TUs with gene (or start of gene)')
    parser.add_argument('-w', '--wig',
        type=int,
        default=1,
        help=f'Index of wig type to use (0-{len(util.WIG_FILE_TEMPLATES)-1}, default: 1)')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

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
    wigs = util.load_wigs(wig_index=args.wig)
    total_reads = np.zeros((2, util.GENOME_SIZE))
    total_reads[0, :] = wigs[0, :] + wigs[2, :]
    total_reads[1, :] = wigs[1, :] + wigs[3, :]
    print('Calculating z scores...')
    z_peak = z_peak_score(wigs, half_width_peak, gap_peak, threshold)
    z_step = z_step_score(total_reads, half_width_step, gap_step, threshold)

    z_step_neg = -z_step
    z_step_pos = z_step
    z_step_neg[z_step_neg < 0] = 0
    z_step_pos[z_step_pos < 0] = 0

    # Load genome info
    genes, locus, all_starts, all_ends = util.load_genome()

    # Identify separate regions of genes for the forward strand
    if args.fwd:
        print('Identifying forward regions...')
        reads = wigs[0, :] + wigs[2, :]
        ma_reads = np.convolve(reads, conv, 'same')
        no_reads = np.where(ma_reads < ma_min_reads)[0]

        starts = all_starts[all_starts > 0]
        ends = all_ends[all_ends > 0]

        region_starts, region_ends = identify_regions(starts, ends, no_reads, gene_pad, ma_pad)
        start_peak = z_peak[2, :]
        end_peak = z_peak[0, :]
        start_step = z_step_neg[0, :]
        end_step = z_step_pos[0, :]
        start_comb = start_peak * start_step
        end_comb = end_peak * end_step
        tu_starts, tu_ends = identify_tu_locations(region_starts, region_ends, start_comb, end_comb)
        tu_genes = identify_tu_genes(tu_starts, tu_ends, genes, all_starts, all_ends)

        # Find regions and TUs for a specified gene or set of genes
        if args.gene:
            for i, tus in enumerate(tu_genes):
                skip = False
                for tu in tus:
                    for g in tu.split(':'):
                        if g.startswith(args.gene):
                            print(i, tus)
                            skip = True
                            break

                    if skip:
                        break

        # Plot forward strand regions
        if args.plot:
            print('Plotting forward regions...')
            labels = ['z peak (start)', 'z peak (end)', 'z step (start)', 'z step (end)', 'combined (start)', 'combined (end)']
            for i, (start, end) in enumerate(zip(region_starts, region_ends)):
                util.plot_tus(start, end, genes, all_starts, all_ends, wigs, tu_genes[i],
                    path=os.path.join(util.OUTPUT_DIR, f'fwd_{i}_tus{args.label}.png'))
                if args.plot_z:
                    scores = np.vstack((
                        start_peak[start:end], end_peak[start:end],
                        start_step[start:end], end_step[start:end],
                        start_comb[start:end], end_comb[start:end],
                        ))
                    util.plot_reads(start, end, genes, all_starts, all_ends, wigs, scores=scores,
                        score_labels=labels, path=os.path.join(util.OUTPUT_DIR, f'fwd_{i}{args.label}.png'))

    # Identify separate regions of genes for the reverse strand
    if args.rev:
        print('Identifying reverse regions...')
        reads = wigs[1, :] + wigs[3, :]
        ma_reads = np.convolve(reads[::-1], conv, 'same')
        no_reads = np.where(ma_reads < ma_min_reads)[0]

        starts = -all_starts[all_starts < 0][::-1]
        ends = -all_ends[all_ends < 0][::-1]

        region_starts, region_ends = identify_regions(starts, ends, no_reads, gene_pad, ma_pad)
        start_peak = z_peak[3, ::-1]
        end_peak = z_peak[1, ::-1]
        start_step = z_step_pos[1, ::-1]
        end_step = z_step_neg[1, ::-1]
        start_comb = start_peak * start_step
        end_comb = end_peak * end_step
        tu_starts, tu_ends = identify_tu_locations(region_starts, region_ends, start_comb[::-1], end_comb[::-1], rev=True)
        tu_genes = identify_tu_genes(tu_starts, tu_ends, genes, all_starts, all_ends)[::-1]

        # Find regions and TUs for a specified gene or set of genes
        if args.gene:
            for i, tus in enumerate(tu_genes):
                skip = False
                for tu in tus:
                    for g in tu.split(':'):
                        if g.startswith(args.gene):
                            print(i, tus)
                            skip = True
                            break

                    if skip:
                        break

        # Plot reverse strand regions
        if args.plot:
            print('Plotting reverse regions...')
            labels = ['z peak (start)', 'z peak (end)', 'z step (start)', 'z step (end)', 'combined (start)', 'combined (end)']
            for i, (start, end) in enumerate(zip(region_starts[::-1], region_ends[::-1])):
                util.plot_tus(-start, -end, genes, all_starts, all_ends, wigs, tu_genes[i],
                    path=os.path.join(util.OUTPUT_DIR, f'rev_{i}_tus{args.label}.png'))
                if args.plot_z:
                    scores = np.vstack((
                        start_peak[-end:-start], end_peak[-end:-start],
                        start_step[-end:-start], end_step[-end:-start],
                        start_comb[-end:-start], end_comb[-end:-start],
                        ))
                    util.plot_reads(-start, -end, genes, all_starts, all_ends, wigs, scores=scores,
                        score_labels=labels, path=os.path.join(util.OUTPUT_DIR, f'rev_{i}{args.label}.png'))
