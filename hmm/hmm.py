
from __future__ import division

import numpy as np
import os
import csv
from matplotlib import pyplot as plt

WIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wigs_U000962')
WIG_FILE = os.path.join(WIG_DIR, 'GSM2971252_Escherichia_coli_WT_Rend_seq_5_exo_MOPS_comp_25s_frag_pooled_{}_no_shadow.wig')
WIG_STRANDS = ['3f', '3r', '5f', '5r']

GENOME_SIZE = 4639675

def load_wigs(wig_file, strands, genome_size):
    '''
    Loads read data from .wig files stored in WIG_DIR

    Inputs:
        wig_file (str) - template for the wig files to be read (gets formatted for each strand)
        strands (list of str) - each strand name that gets inserted into the wig_file string
            (order determines order of first dimension of returned array)
        genome_size (int) - base pair length of the genome

    Returns 2D np.array of floats (4 x genome_size) for reads on each strand
    '''

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
                    print exc

            index = np.array(d[:, 0], dtype=int) - 1
            reads = np.array(d[:, 1], dtype=float)

            data[i, index] = reads

    return data

if __name__ == '__main__':
    reads = load_wigs(WIG_FILE, WIG_STRANDS, GENOME_SIZE)

    # f1
    start = 1
    end = 5233
    sigma_on = 0.5
    sigma_off = 0.5

    # f1
    start = 16547
    end = 21180
    sigma_on = 0.6
    sigma_off = 1.

    # # f33 and f34
    # start = 164524
    # # end = 167422
    # # start = 167254
    # end = 175106
    # sigma_on = 0.6
    # sigma_off = 1.

    # # f39 and f40
    # start = 183610
    # # end = 194775
    # # start = 194733
    # end = 214290
    # sigma_on = 0.5
    # sigma_off = 0.5

    # # f823
    # start = 3427802
    # end = 3436017
    # sigma_on = 0.6
    # sigma_off = 1.

    total_reads = reads[0, start:end] + reads[2, start:end]

    p = 3000 / GENOME_SIZE
    pi = np.array([0.5, 0.5])

    A = np.array([[1-p, p], [p, 1-p]])
    Bon = np.array([[np.exp(-1/(2*sigma_on**2)), 0], [0, 1]]) / np.sqrt(2*np.pi*sigma_on**2)
    Boff = np.array([[1, 0], [0, np.exp(-1/(2*sigma_off**2))]]) / np.sqrt(2*np.pi*sigma_off**2)

    Sf = np.zeros(end-start)
    alpha = np.zeros((2, end-start))
    alpha[:, 0] = pi
    for i, read in enumerate(total_reads[1:]):
        if read > 0:
            B = Bon
        else:
            B = Boff

        probs = np.dot(np.dot(B, A.T), alpha[:, i].T)

        alpha[:, i+1] = probs / np.sum(probs)
        Sf[i+1] = alpha[1, i+1]

    Sr = np.zeros(end-start)
    beta = np.zeros((2, end-start))
    beta[:, 0] = np.ones_like(pi)
    for i, read in enumerate(total_reads[-1:0:-1]):
        if read > 0:
            B = Bon
        else:
            B = Boff

        probs = np.dot(np.dot(beta[:, i].T, A), B)

        beta[:, i+1] = probs / np.sum(probs)
        Sr[i+1] = beta[1, i+1]

    smoothed = alpha * beta[:, ::-1]
    gamma = smoothed / np.sum(smoothed, axis=0)

    plt.figure(figsize=(20,10))
    plt.plot(range(start,end), np.log10(total_reads))
    plt.plot(range(start,end), Sf, 'r')
    plt.plot(range(start,end), Sr[::-1], 'g')
    plt.plot(range(start,end), gamma[1, :] > 0.1, 'k')
    plt.show()
