'''
Script to convert .fasta file from http://rocaplab.ocean.washington.edu/tools/convert-genbank-to-fasta/
to the .faa format expected for the isoforms matlab script.

Expects header lines to be separated by | with values locus_tag, gene, protein, location
'''

import re


INPUT_FILE = 'U00096.2.fasta'
OUTPUT_FILE = 'U00096.2.faa'

HEADERS = ['locus_tag', 'gene', 'protein', 'location']

with open(INPUT_FILE, 'r') as f:
	data = list(f)

with open(OUTPUT_FILE, 'w') as f:
	for line in data:
		if line.startswith('>'):
			# ignore > at beginning of line for processing
			values = line[1:].split(' | ')
			f.write('>')

			# write each header value with appropriate formatting
			for value, header in zip(values, HEADERS):
				if header == 'location':
					location = value.split(' ')
					location[0] = re.sub(':', '..', location[0])
					if location[1] == 'Reverse\n':
						value = 'complement({})'.format(location[0])
					else:
						value = location[0]
				f.write('[{}={}] '.format(header, value))
			f.write('\n')
		else:
			f.write(line)
