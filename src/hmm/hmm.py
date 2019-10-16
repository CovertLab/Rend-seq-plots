'''
Forward-backward method based on Krishnamurthy and Chung. "Signal Processing Based on Hidden Markov
Models for Extracting Small Channel Currents". 2006.
'''


from __future__ import division

import json
import numpy as np
import os
import csv
from matplotlib import pyplot as plt

# import ipdb; ipdb.set_trace()
WIG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data/wigs')
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
                    print(exc)

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
    alpha[0, 0] = 1
    for i, read in enumerate(signal[1:]):
        B = np.diag(calc_gaussian(read, levels, variance))

        probs = np.dot(np.dot(B, A), alpha[:, i].T)
        total = np.sum(probs)
        if total == 0:
            alpha[:, i+1] = alpha[:, i]
        else:
            alpha[:, i+1] = probs / total

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

    n_levels = len(levels)
    n_points = len(signal)

    A = A.T
    Bs = np.zeros((n_levels, n_levels, n_points))
    beta = np.zeros((n_levels, n_points))
    beta[:, 0] = 1
    for i, read in enumerate(signal[-1:0:-1]):
        B = np.diag(calc_gaussian(read, levels, variance))

        probs = np.dot(np.dot(beta[:, i].T, A), B)
        total = np.sum(probs)
        if total == 0:
            beta[:, i+1] = beta[:, i]
            Bs[:, :, -(i+1)] = Bs[:, :, -i]  # could be issue if i = 0
        else:
            beta[:, i+1] = probs / total
            Bs[:, :, -(i+1)] = B

    return beta[:, ::-1], Bs

if __name__ == '__main__':
    reads = load_wigs(WIG_FILE, WIG_STRANDS, GENOME_SIZE)

    ## information for specific test cases
    # n_levels should be number of genes + 2 (for level at 0 and for peak level)

    # f1
    start = 1
    end = 5233
    sigma_on = 0.5
    n_levels = 4
    # total_reads = np.array(json.load(open('hmm/f1_reads.json')))
    total_reads = reads[0, start:end] + reads[2, start:end]
    
    ''' with open('f1_reads.json', 'w') as f:
        json.dump(list(total_reads), f)
    # json.dump(list(total_reads), open('hmm/f1_reads.json', 'w'))

    # # f5
    # start = 16547
    # end = 21180
    # sigma_on = 0.6
    # n_levels = 4
    # total_reads = np.array(json.load(open('hmm/f5_reads.json')))
    # # total_reads = reads[0, start:end] + reads[2, start:end]
    # # json.dump(list(total_reads), open('hmm/f5_reads.json', 'w'))

    # # f33 and f34
    # start = 164524
    # # end = 167422
    # # start = 167254
    # end = 175106
    # sigma_on = 0.6
    # n_levels = 7
    # total_reads = np.array(json.load(open('hmm/f33-34_reads.json')))
    # # total_reads = reads[0, start:end] + reads[2, start:end]
    # # json.dump(list(total_reads), open('hmm/f33-34_reads.json', 'w'))

    # # f39 and f40
    # start = 183610
    # # end = 194775
    # # start = 194733
    # end = 214290
    # sigma_on = 0.5
    # n_levels = 22
    # total_reads = np.array(json.load(open('hmm/f39-40_reads.json')))
    # # total_reads = reads[0, start:end] + reads[2, start:end]
    # # json.dump(list(total_reads), open('hmm/f39-40_reads.json', 'w'))

    # # f823
    # start = 3427802
    # end = 3436017
    # sigma_on = 0.6
    # n_levels = 6
    # total_reads = np.array(json.load(open('hmm/f823_reads.json')))
    # # total_reads = reads[0, start:end] + reads[2, start:end]
    # # json.dump(list(total_reads), open('hmm/f823_reads.json', 'w'))

    ## two level HMM test case with forward-backward algorithm
    # total_reads = reads[0, start:end] + reads[2, start:end]
    #
    # pi = np.array([0.5, 0.5])
    # p = 2500 / GENOME_SIZE
    # A = np.array([[1-p, p], [p, 1-p]])
    #
    # alpha = forward_algo([1 if x > 0 else 0 for x in total_reads], pi, A, np.array([0, 1]), sigma**2)
    # beta, _ = backward_algo([1 if x > 0 else 0 for x in total_reads], A, np.array([0,1]), sigma**2)
    #
    # smoothed = alpha * beta
    # gamma = smoothed / np.sum(smoothed, axis=0)
    #
    # plt.figure(figsize=(20,10))
    # plt.plot(range(start,end), np.log10(total_reads))
    # plt.plot(range(start,end), alpha[1, :], 'r')
    # plt.plot(range(start,end), beta[1, :], 'g')
    # plt.plot(range(start,end), gamma[1, :] > 0.1, 'k')
    # plt.show()

    ## EM method for determining HMM levels
    sigma = np.std(total_reads)
    levels = sigma * np.ones(n_levels)
    levels[0] = 0
    levels[-1] = 0 #np.array([0, 6, 6, 0])
    # sigma = np.sqrt(np.max(total_reads))
    # levels = np.linspace(0, np.max(total_reads), n_levels)
    print('Sigma: {}\nLevels: {}'.format(sigma, levels))

    n_iterations = 20
    pi = np.ones(n_levels) / n_levels
    p = 2500 / GENOME_SIZE
    # transition matrix assumes 1 transition to and from a level starting and ending at 0
    A = (1-p) * np.eye(n_levels)
    A[:-1, 1:] += p / 2 * np.eye(n_levels - 1)

    for i in range(n_iterations):
        alpha = forward_algo(total_reads, pi, A, levels, sigma**2)
        beta, B = backward_algo(total_reads, A, levels, sigma**2)

        alpha_diag = np.array([np.diag(col) for col in alpha.T])
        beta_diag = np.array([np.diag(col) for col in beta.T])

        a_gamma = np.dot(alpha_diag[:-1,:,:], A)
        b_gamma = beta_diag[1:,:,:] * B.T[1:,:,:]
        smoothed = np.array([np.dot(a, b) for a, b in zip(a_gamma, b_gamma)])
        sums = np.sum(np.sum(smoothed, axis=1), axis=1).reshape(-1, 1, 1)
        gamma = smoothed / sums

        gamma_ik = np.sum(gamma, axis=2)
        gamma_ij = np.sum(gamma, axis=0)
        gamma_i = np.sum(gamma_ij, axis=1)

        # calculate updated A matrix (eq 19.40)
        numerator = gamma_ij
        denominator = np.zeros_like(numerator)
        for i in range(len(denominator)):
            denominator[i, :] = gamma_i[i]
        A_new = numerator / denominator

        # calculate updated state levels (eq 19.41)
        numerator = np.sum(gamma_ik * total_reads[1:].reshape(-1, 1), axis=0)
        denominator = np.sum(gamma_ik, axis=0)
        q_new = numerator / denominator

        # calculate updated standard deviation values (eq 19.42)
        sigma_new = np.sqrt(1 / len(total_reads) * np.sum(gamma_ik * (np.tile(total_reads[1:], (len(q_new),1)).T - q_new)**2))

        A = A_new
        levels = q_new
        levels[0] = 0  # know level must be at 0 for no reads, doesn't impact convergence or final values significantly
        levels[-1] = 0
        sigma = sigma_new

    print(levels)
    assigned_levels = np.array([levels[np.argmax(np.diag(x))] for x in gamma])

    # Potentially useful to use moving average to calculate levels
    window = 10
    clipped_index = int((window-1) / 2)  # For proper indexing since data is lost
    ma = np.convolve(total_reads, np.ones((window,))/window, 'valid')'''
    
    x = range(start,end)
    g1 = reads[0, start:end]
    g2 = reads[2, start:end]
    fig, ax = plt.subplots()
    plt.scatter(g1,g2)
    plt.title('Scatter plot for 3f vs 5f')
    # plt.legend(loc='4')
    plt.xlabel('3f')
    plt.ylabel('5f')

    ax.plot(ax.get_xlim(), ax.get_ylim())    
    plt.savefig("scatterplot-35.pdf")
    # plt.savefig("scatterplot.pdf")

    '''ax.hist(g1, color='blue')
    ax.hist(g2, color='orange')
    plt.show()
    # plt.savefig('Histogram')'''
    import ipdb; ipdb.set_trace()


'''
    # Plot histrogram of level distributions
    plt.figure()
    for i, level in enumerate(levels):
        # Filter high outliers
        reads = total_reads[1:][assigned_levels == level]
        mean = np.mean(reads)
        std = np.std(reads)
        filtered_reads = reads[reads < mean + 5*std]
        n_bins = int(np.max(filtered_reads) - np.min(filtered_reads))

        # # Plot distributions within a given level
        # if i == 1:
        #     size = 20
        #     plt.figure(figsize=(5,10))
        #     n_plots = int(len(filtered_reads) / size)
        #     for j in range(n_plots):
        #         ax = plt.subplot(n_plots, 1, j+1)
        #         r = filtered_reads[j*size:(j+1)*size]
        #         ax.hist(r, bins=range(50))
        #     plt.show()

        ax = plt.subplot(n_levels, 1, i+1)
        ax.hist(filtered_reads, bins=n_bins)
        mask = assigned_levels[clipped_index:-clipped_index] == level
        ax.hist(ma[mask][ma[mask] < mean + 3*std], color='r', alpha=0.5)
    plt.show()

    # Plot sample value vs moving average value for different levels
    # Could be useful for Gaussian Dicriminant Analysis to classify levels
    plt.figure(figsize=(5,10))
    for i, level in enumerate(levels):
        mask = assigned_levels[clipped_index:-clipped_index] == level
        ax = plt.subplot(n_levels, 1, i+1)
        ax.plot(ma[~mask], total_reads[clipped_index+1:-clipped_index][~mask], 'o', color='gray', alpha=0.1)
        ax.plot(ma[mask], total_reads[clipped_index+1:-clipped_index][mask], 'rx', alpha=0.3)

        # # Focus on more populated area for f1 data
        # ax.set_xlim([0, 15])
        # ax.set_ylim([0, 15])
    plt.show()

    # Plot levels and read data
    plt.figure(figsize=(20,10))
    plt.plot(range(start,end), np.log10(total_reads+1))
    plt.plot(range(start,end)[1:], np.log10(assigned_levels+1))
    plt.show()
    '''
