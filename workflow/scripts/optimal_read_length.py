#!/usr/bin/env python3
from __future__ import print_function
import sys, glob, json, re


def find_optimal_read_length(rl, dbrl):
	# Get the best read length given a list of available STAR Indices
	return next(x[1] for x in enumerate(dbrl) if x[1] >= rl)


if __name__ == '__main__':
	# Parse required positional args
	try:
		readlength = sys.argv[1] # Max read length of all samples calculated from FastQC
		stardir = sys.argv[2]    # PATH to STAR Indices for different read lengths
	except IndexError:
		print('Failed to provide all required positional args!')
		print('Example Usage: python {} QC/readlength.txt /data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_30/STAR/2.7.0f/genes-'.format(sys.argv[0]))
		sys.exit(1)

	# Get max read length of sample
	my_read_length=int(open(readlength).readlines()[0].strip())-1

	# Find all STAR Indice read lengths
	star_read_lengths=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]), glob.glob(stardir+'*/'))))

	myrl = find_optimal_read_length(rl=my_read_length, dbrl=star_read_lengths)

	print(myrl)
