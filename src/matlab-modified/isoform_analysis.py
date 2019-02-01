'''
Analysis script to further process output from Rend_seq_isoform_header_20180130.m

- Removes isoforms that cover the same gene or set of genes (with different start/end points)
- Combines the level for duplicate isoforms
'''

import numpy as np
import csv

ISOFORM_FILE = 'isoforms.tsv'
HEADERS = [
	'Operon',
	'# genes',
	'Strand',
	'Start',
	'Stop',
	'Level',
	'Genes',
	]


# Load data from output from Rend_seq_isoforms_header.m
with open(ISOFORM_FILE, 'r') as f:
	reader = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONNUMERIC, quotechar='"')

	data = list(reader)

# Remove headers from data
headers = data[0]
data = data[1:]
operon_idx = headers.index('Operon')
n_genes_idx = headers.index('# Genes')
genes_idx = headers.index('Genes')
level_idx = headers.index('Level')
strand_idx = headers.index('Strand')
frame_idx = headers.index('Frame')

# Data structures for data of interest
n_genes_per_operon = {}
tus_in_operon = {}
level_per_tu = {}
operon_frame = {}

# Process data
for row in data:
	operon_num = int(row[operon_idx])
	n_genes = int(row[n_genes_idx])
	genes = row[genes_idx]
	level = row[level_idx]
	strand = row[strand_idx]
	frame = int(row[frame_idx])

	n_genes_per_operon[operon_num] = n_genes

	gene_set = tus_in_operon.get(operon_num, set())
	if genes in gene_set:
		level_per_tu[genes] += level
	else:
		level_per_tu[genes] = level

		gene_set.add(genes)
		tus_in_operon[operon_num] = gene_set

	operon_frame[operon_num] = '{}{}'.format(strand, frame)

# Check for underdetermined operons
for operon in n_genes_per_operon:
	n_tus = len(tus_in_operon[operon])
	n_genes = n_genes_per_operon[operon]
	if n_tus > n_genes:
		print('Underdetermined operon {} ({}): {} genes in {} TUs'.format(operon, operon_frame[operon], n_genes, n_tus))


import ipdb; ipdb.set_trace()
# print('{} total operons'.format(data[-1, 0]))
# print('{} operons contain a single gene'.format(sum(data[:, 1] == 1)))
