#! /usr/bin/env python
'''
Identifies isoforms based on data in .wig files.
'''

import os

import numpy as np

import util


label = '_0.05'


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
    threshold = 0.05
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
    print('Calculating z scores...')
    z_peak, z_step = util.z_score_statistics(wigs, half_width_peak=half_width_peak,
        half_width_step=half_width_step, threshold=threshold, gap_peak=gap_peak, gap_step=gap_step)

    z_peak[z_peak < 1] = 1
    z_step[z_step < 1] = 1

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
    for i, (start, end) in enumerate(zip(real_starts, real_ends)):
        util.plot_reads(start, end, genes, all_starts, all_ends, wigs,
            path=os.path.join(util.OUTPUT_DIR, f'fwd_{i}.png'))
        util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=start_peak[start:end],
            path=os.path.join(util.OUTPUT_DIR, f'fwd_{i}_starts_peak{label}.png'))
        util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=end_peak[start:end],
            path=os.path.join(util.OUTPUT_DIR, f'fwd_{i}_ends_peak{label}.png'))

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
    for i, (start, end) in enumerate(zip(real_starts[::-1], real_ends[::-1])):
        util.plot_reads(-start, -end, genes, all_starts, all_ends, wigs,
            path=os.path.join(util.OUTPUT_DIR, f'rev_{i}.png'))
        util.plot_reads(-start, -end, genes, all_starts, all_ends, wigs, fit=start_peak[-end:-start],
            path=os.path.join(util.OUTPUT_DIR, f'rev_{i}_starts_peak{label}.png'))
        util.plot_reads(-start, -end, genes, all_starts, all_ends, wigs, fit=end_peak[-end:-start],
            path=os.path.join(util.OUTPUT_DIR, f'rev_{i}_ends_peak{label}.png'))

    # Explore data for an example region
    # combined = z_peak * z_step
    # combined1 = z_peak[0, :] * z_step[2, :]
    # combined2 = z_peak[2, :] * z_step[0, :]
    # z_peak[z_peak < 1] = 1
    # z_step[z_step < 1] = 1
    # combined[combined < 1] = 1
    # combined1[combined1 < 1] = 1
    # combined2[combined2 < 1] = 1
    # start = 1
    # end = 6000
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=z_peak[0, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_peak0.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=z_peak[2, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_peak2.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=z_step[0, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_step0.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=z_step[2, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_step2.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=combined[0, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_z0.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=combined[2, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_z2.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=combined1[start:end], path=os.path.join(util.OUTPUT_DIR, '{}_combined1.png'.format(name)))
    # util.plot_reads(start, end, genes, all_starts, all_ends, wigs, fit=combined2[start:end], path=os.path.join(util.OUTPUT_DIR, '{}_combined2.png'.format(name)))
    import ipdb; ipdb.set_trace()
