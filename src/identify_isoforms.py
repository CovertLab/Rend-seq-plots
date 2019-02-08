#! /usr/bin/env python
'''
Identifies isoforms based on data in .wig files.
'''

import os

import util


if __name__ == '__main__':
    half_width_peak = 100
    half_width_step = 100
    threshold = 0.25
    gap_peak = 2  # not used unless filter is changed in z_score_statistics
    gap_step = 3

    wigs = util.load_wigs()
    z_peak, z_step = util.z_score_statistics(wigs, half_width_peak=half_width_peak,
        half_width_step=half_width_step, threshold=threshold, gap_peak=gap_peak, gap_step=gap_step)

    # Explore data for an example region
    z_peak[z_peak < 1] = 1
    z_step[z_step < 1] = 1
    combined = z_peak * z_step
    combined1 = z_peak[0, :] * z_step[2, :]
    combined2 = z_peak[2, :] * z_step[0, :]
    genes, locus, starts, ends = util.load_genome()
    start = 1
    end = 6000
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_peak[0, start:end], path=os.path.join(util.OUTPUT_DIR, 'peak0.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_peak[2, start:end], path=os.path.join(util.OUTPUT_DIR, 'peak2.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_step[0, start:end], path=os.path.join(util.OUTPUT_DIR, 'step0.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_step[2, start:end], path=os.path.join(util.OUTPUT_DIR, 'step2.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined[0, start:end], path=os.path.join(util.OUTPUT_DIR, 'z0.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined[2, start:end], path=os.path.join(util.OUTPUT_DIR, 'z2.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined1[start:end], path=os.path.join(util.OUTPUT_DIR, 'combined1.png'))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined2[start:end], path=os.path.join(util.OUTPUT_DIR, 'combined2.png'))
    import ipdb; ipdb.set_trace()
