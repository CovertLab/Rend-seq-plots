#! /usr/bin/env python
'''
Identifies isoforms based on data in .wig files.
'''

import os

import util

name = 'test'

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
    combined = z_peak * z_step
    combined1 = z_peak[0, :] * z_step[2, :]
    combined2 = z_peak[2, :] * z_step[0, :]
    z_peak[z_peak < 1] = 1
    z_step[z_step < 1] = 1
    combined[combined < 1] = 1
    combined1[combined1 < 1] = 1
    combined2[combined2 < 1] = 1
    genes, locus, starts, ends = util.load_genome()
    start = 1
    end = 6000
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_peak[0, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_peak0.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_peak[2, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_peak2.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_step[0, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_step0.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=z_step[2, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_step2.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined[0, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_z0.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined[2, start:end], path=os.path.join(util.OUTPUT_DIR, '{}_z2.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined1[start:end], path=os.path.join(util.OUTPUT_DIR, '{}_combined1.png'.format(name)))
    util.plot_reads(start, end, genes, starts, ends, wigs, fit=combined2[start:end], path=os.path.join(util.OUTPUT_DIR, '{}_combined2.png'.format(name)))
    import ipdb; ipdb.set_trace()
