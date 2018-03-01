#!/usr/bin/python
"""
Usage:
python trim_adapter.py 1 read_1.fasta.gz read_1_trimmed.fasta.gz
"""

import gzip
import sys

is_left = True if sys.argv[1] == '1' else False
filepath_in = sys.argv[2]
filepath_out = sys.argv[3]

with gzip.open(filepath_in, 'rb') as f:
	data = f.read().split("\n")

total_count = len(data)/4

with gzip.open(filepath_out, 'wb') as f:
	for i in range(total_count):
		header = data[i*4]
		read = data[i*4+1]
		plus = data[i*4+2]
		quality = data[i*4+3]
		
		if is_left:
			adapter = header.split(':')[-1].split('+')[0]
		else:
			adapter = header.split(':')[-1].split('+')[1].strip('N')
		trim_indx = read.find(adapter)
		
		f.write('%s\n%s\n%s\n%s\n' % (header, read[:trim_indx], plus, quality[:trim_indx]))
