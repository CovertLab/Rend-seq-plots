
from __future__ import division

import numpy as np
import os
import csv
from matplotlib import pyplot as plt

WIG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
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

# Functions for HMM algorithm
def calc_gaussian(signal, levels, variance):
    '''
    Calculates the probability based on a normal distribution for b_i(y_k)
    Based on Krishnamurthy eq 19.13

    Inputs:
        signal (float) - y_k, level of signal at position k
        levels (array of floats) - q, level of signal for each possible state
        variance (float) - sigma_w^2, variance for the signal

    Returns array of floats for the probability of each possible state (does not need to sum to 1)
    '''

    return np.exp(-0.5*(signal - levels)**2 / variance) / np.sqrt(2*np.pi*variance)

def forward_algo(signal, pi, A, levels, variance):
    '''
    Calculates alpha for the forward algorithm
    Calculates B based on the levels and variance
    Based on Krishnamurthy eq 19.27

    Inputs:
        signal (array of floats) - y, signal at each position k
        pi (array of floats) - initial value for alpha
        A (2D array of floats) - transition probability matrix
        levels (array of floats) - q, level of signal for each possible state
        variance (float) - sigma_w^2, variance for the signal

    Returns 2D array of floats for alpha (n states by k positions)
    '''

    A = A.T
    alpha = np.zeros((len(levels), len(signal)))
    alpha[:, 0] = pi
    for i, read in enumerate(signal[1:]):
        B = np.diag(calc_gaussian(read, levels, variance))

        probs = np.dot(np.dot(B, A), alpha[:, i].T)
        alpha[:, i+1] = probs / np.sum(probs)

    return alpha

def backward_algo(signal, A, levels, variance):
    '''
    Calculates beta for the backward algorithm
    Calculates B based on the levels and variance
    Initial values of beta are set to 1
    Based on Krishnamurthy eq 19.35

    Inputs:
        signal (array of floats) - y, signal at each position k
        A (2D array of floats) - transition probability matrix
        levels (array of floats) - q, level of signal for each possible state
        variance (float) - sigma_w^2, variance for the signal

    Returns 2D array of floats for beta (n states by k positions)
    '''

    A = A.T
    beta = np.zeros((len(levels), len(signal)))
    beta[:, 0] = 1
    for i, read in enumerate(signal[-1:0:-1]):
        B = np.diag(calc_gaussian(read, levels, variance))

        probs = np.dot(np.dot(beta[:, i].T, A), B)
        beta[:, i+1] = probs / np.sum(probs)

    return beta

if __name__ == '__main__':
    reads = load_wigs(WIG_FILE, WIG_STRANDS, GENOME_SIZE)

    # information for specific test cases
    # f1
    start = 1
    end = 5233
    sigma_on = 0.5
    sigma_off = 0.5

    # f5
    # start = 16547
    # end = 21180
    # sigma_on = 0.6
    # sigma_off = 1.

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

    # two level HMM test case with forward-backward algorithm
    total_reads = reads[0, start:end] + reads[2, start:end]

    pi = np.array([0.5, 0.5])
    p = 3000 / GENOME_SIZE
    A = np.array([[1-p, p], [p, 1-p]])

    alpha = forward_algo([1 if x > 0 else 0 for x in total_reads], pi, A, np.array([0, 1]), sigma_on**2)
    beta = backward_algo([1 if x > 0 else 0 for x in total_reads], A, np.array([0,1]), sigma_on**2)

    smoothed = alpha * beta[:, ::-1]
    gamma = smoothed / np.sum(smoothed, axis=0)

    plt.figure(figsize=(20,10))
    plt.plot(range(start,end), np.log10(total_reads))
    plt.plot(range(start,end), alpha[1, :], 'r')
    plt.plot(range(start,end), beta[1, ::-1], 'g')
    plt.plot(range(start,end), gamma[1, :] > 0.1, 'k')
    plt.show()
