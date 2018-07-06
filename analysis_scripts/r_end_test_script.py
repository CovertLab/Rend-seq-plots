# Import Packages
import glob
import numpy as np
#from numpy.fft import fft, ifft, fftshift
#import matplotlib.pyplot as plt
import os
import pandas as pd
import re
#import scipy as sp

import pyfftw.interfaces.scipy_fftpack as fftw

# Make these variables global?
data_subfolder = 'wigs_U000962'
data_file_name = "GSM2971252_Escherichia_coli_WT_Rend_seq_5_exo_MOPS_comp_25s_frag_pooled_*_no_shadow.wig"
genome_annotation_file_name = "U00096.2.faa"
genome_size = 4639675 # length of genome in base pairs

# Tunable parameters:

half_width_z_score = 100   # half width of window for averaging for peak z score
half_width_step = 100.      # half width of window for averaging for step z score
average_threshold = 0.25   # average read density (read/nt) threshold for consideration
#gap_z = 2.                  # gap left out (both sides) of central position for peak z score 
gap_step = 3.   

# Add something else to make sure fft is calculated on an even array

#Functions for data loading and parsing
def load_wig_data(subfolder, data_file_name, genome_len):
    """
    The purpose of this function is to go create a data file for all the 
    wig sequencing data. The data is is in the form of a text file (.wig) and contains position 
    (first column) and count (second column) data.
    
    This function will simply take in the data files, and put them into an array for use by downstream functions.
    
    Requires:

    -Subfolder where files are located.
    -Name of data file containing a wildcard in the place of the directional info for that sequencing run.
    
    Returns:
    
    - dense_data_array: a 3 dimensional array.
        First dimension: Based on seq file.
        Second dimension: counts per position
        All positions through the lenght of the genome are recorded.
    - file_paths: a list containing all of the file paths 
    - seq_directions: a list containing the sequencing direction of each file,
        in the order it was pulled to create all_data. The directions are either 3' or 
        5' in either forward or reverse direction.
    """
    #parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    data_directory = os.path.join(os.getcwd(), 'data', subfolder, data_file_name)
    file_paths = []
    dense_data_array = np.zeros((4, genome_len + 1, 1))
    count = 0
    for file_path in glob.glob(data_directory):
        file_paths.append(file_path)
        with open (file_path, 'rt') as f:
            for line in f:
                if line[0].isdigit():
                    x, y = line.split()
                    dense_data_array[count, int(x) - 1] = float(y)
        count += 1
    seq_directions = get_seq_direction(file_paths)
    return file_paths, dense_data_array, seq_directions


def get_seq_direction(file_names):
    """
    Purpose:
    - To extract from the file names the seq direction of the imported wig file. 
        This is to avoid having to hard-code the directionality in the order the files 
        are pulled, although the direction options do need to be hard-coded, and 
        need to be present in the file name.
    Requires:
    - List of all files that need to be evaluated (this list is output in the function
        load_wig_data)
    - List of all the possible directions to be evaluated
        This list of directions should be the same regardless of the 
        sample. In this function these directions are hard coded. See possible_directions.
    Returns:
    -Read direction as a list of strings based on the order the files were imported 
        and are consequently stored.
    - The positions are listed in the order that the files were pulled and extracted
    in the function load_wig_data(). So the positions in this list reflect the 
    order of that data.
    """
    possible_directions = ['5f', '5r', '3f', '3r']
    directions = []
    for name in file_names:
        for direction in possible_directions:
            if direction in name.lower():
                directions.append(direction)
    return directions

def simplify_genome_annotation(annotation_file_name):
    """
    The purpose of this function is to pull out genome annotation information from the 
    original .faa file.

    Requires:
        -.faa genome annotation file in /genome/ subdirectory.

    Returns:
        -List where for each gene the the locus tag, gene name, gene name,
            gene direction, start position, and end position are recorded
            For example:
            ['b0001', 'thrL', '+', 190.0, 255.0]

    Note:
        -The annotation file has to be annotated in a very specific format.
        For example:
        >[locus_tag=b0585] [gene=fes] [protein=enterobactin/ferric enterobactin esterase] [location=612038..613162]
        -If file is not in correct format, update code below, or restructure file.
    """
    #parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    file_path = os.path.join(os.getcwd(), 'genome', genome_annotation_file_name)
    line_marker = '>'
    annotation_list = []
    with open (file_path, 'rt') as in_file:  
        for line in in_file:
            if line_marker in line:
                annotation_list.append(parse_single_gene(line))
    return annotation_list

def parse_single_gene(single_line):
    """
    Single text line (corresponding to a single gene) originally extracted from gene 
    anotation file in the function simplify_genome_annotation(), and parses it to 
    give extract the anotation information.
    Requires:
        - Single gene text line extracted from the genome annotation file.
    
    Returns:
        -A list contining the locus tag, gene name, gene direction ('+' or '-')(as strings), 
            and gene start and stop for each gene.
    Note:
        - For get_gene_direction_info the data supplied are taken from
            gene_info_not_parsed[-1], since sometimes protein descriptions
            contain [] and will throw off indexing if choosing the third 
            position. But the location information is always stored last.
    """
    gene_info_not_parsed = splice_string_to_attributes(single_line)
    locus_tag = get_locus_tag(gene_info_not_parsed[0])
    gene_name = get_gene_name(gene_info_not_parsed[1])
    gene_direction, gene_start, gene_stop = get_gene_direction_info(gene_info_not_parsed[-1])
    return [locus_tag, gene_name, gene_direction, gene_start, gene_stop]  

def splice_string_to_attributes(single_line):
    """
    Requires:
        -A line of text from the genome annotation file, 
    Returns:
        - list containing the parts necessary for further down stream string
        splicing. 
        -At this point it simply separates information based on whether it is surrounded by square
        brackets.
    """
    split_string = re.split(r'[\[\]]', single_line)
    strip_list = [gene_data for gene_data in split_string if gene_data.strip()][1:]
    return strip_list

def get_locus_tag(locus_info):
    """
    Requires:
        -String in the form: locus_tag=b0586,
    Returns:
        -Locus tag (string)
    """
    locus = locus_info.split('=')[1]
    return locus

def get_gene_name(gene_info):
    """
    Requires:
        -String in the form: gene=ybdZ,
    Returns:
        -Gene name (string)
    Note:
    Requires that no more than a single gene name is assigned per gene. 
    """
    name = gene_info.split('=')[1]
    return name

def get_gene_direction_info(direction_info):
    """
    Requires:
        - Takes in a string in the form: location=417113..418408 or
        location=complement(414974..416176).
    Returns:
        -If is the complemnent returns the direction as '-', or 
        else returns '+'.
        -Returns the start and stop positions as floats.
    """
    start = [] 
    stop = []
    location_str = direction_info.split('=')[1]
    direction  = ['-' if '(' in location_str else '+'][0]
    
    if direction == "+":
        start, stop = location_str.split('..')
    else:
        indices = [location_str.find(i) for i in ['(', ')']]
        start, stop = location_str[indices[0]+1:indices[1]].split('..')
    return direction , float(start), float(stop)

#Functions for computations


def create_fft_filter(direction, genome_len, half_width):
    filt = np.zeros((genome_len + 1, 1))

    if 'f' in direction:
        filt[genome_len / 2 + 1 : genome_len / 2 + 2 * half_width_z_score] = 1. / (2. * half_width_z_score)
    elif 'r' in direction:
        filt[genome_len / 2 - 2 * half_width_z_score : genome_len / 2 - 1] = 1. / (2. * half_width_z_score)
    else:
        'Fourier tranform filter could not be computed bc the direction of the seq file has not been recorded properly'
    return np.squeeze(filt)

#def create_z_step_filter(data_length, half_width, gap_width):
#    return

def calculate_single_peak_z_score(data_single_series, seq_direct, genome_len, half_width):
    """
    Assumes the single data series has already been squeezed.

    """

    len_data = len(data_single_series)

    print('Making filter')
    filter_for_fft = create_fft_filter(seq_direct, genome_size, half_width)
    print('Performing fft on data')
    fft_data = fft(data_single_series, len_data)
    print('Performing fft on filter')
    fft_filter = fft(filter_for_fft, len_data)
    print('Performing fft data squared')
    fft_data_squared = fft(np.square(data_single_series), len_data)
    print('Getting average data... this will take a bit, be patient')
    average_data = fftshift(ifft((fft_filter.T * fft_data), len_data))
    print('Getting std just keep it together now.')
    std_data = np.sqrt(fftshift(ifft(fft_filter.T * fft_data_squared, len_data)) - np.square(average_data))
    print('Getting z peak scores.... hang in there')
    z_peak = (average_data > average_threshold) * (data_single_series - average_data) / (std_data + (average_data == 0))

    return average_data, z_peak

#def caluculate_single_peak_z_score(data_single_series, half_width, avg_thres, gap_step):
#    return

def calculate_all_z_statistics(all_data, seq_directions, genome_len, half_width_z_score, half_width_step, avg_thresh, gap_step):

    # initialize empty vectors to store data
    all_z_peak = np.zeros((genome_len + 1, 4))
    all_average_data = np.zeros((genome_len + 1, 4))
    all_z_step = np.zeros((genome_len + 1, 4))

    count = 0

    print('Currently on series ' + str(count))

    for count in range(0, len(all_data)):

        [all_average_data[:, count], all_z_peak[:, count]] = calculate_single_peak_z_score(np.squeeze(all_data[count]), 
            seq_directions[count], genome_len, half_width_z_score)
        print('Okay dont freak out... gotta get the z step score. ')
        #all_z_step[:, count] = calculate_single_step_z_score(np.squeeze(all_data[count]), half_width_step, avg_thresh, gap_step)

    return all_z_peak, all_average_data, all_z_step

all_file_names, all_data, seq_directions = load_wig_data(data_subfolder, data_file_name, genome_size)
genome_annotation = simplify_genome_annotation(genome_annotation_file_name)
#import pdb; pdb.set_trace()


#all_computed_statistics = calculate_all_z_statistics(all_data, seq_directions, genome_size, half_width_z_score, half_width_step, average_threshold, gap_step)

single_data_set = np.squeeze(all_data[3])
len_data = len(single_data_set)
direction = seq_directions[3]

filter_for_fft = create_fft_filter(direction, genome_size, half_width_z_score)

fft_data = fftw.fft(single_data_set, len_data)



#ft_filter = np.fft.fft(filter_for_fft, len_data)


#fft_data_squared = np.fft.fft(np.square(single_data_set), len_data) # This step takes ~20 seconds

#average_data = np.fft.fftshift(np.fft.ifft((ft_filter.T * fft_data), len_data))

#std_data = np.sqrt(np.fft.fftshift(np.fft.ifft(ft_filter.T * fft_data_squared, len_data)) - np.square(average_data))

#z_peak = (average_data > average_threshold) * (single_data_set - average_data) / (std_data + (average_data == 0))

