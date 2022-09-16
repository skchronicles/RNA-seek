#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import sys, pandas, os


def create(file, tin_dict, key_index=0, parse_index=4):
	"""Populates the TIN nested dictionary
	@param file <str>: Path to RSEQC output file with TIN values to extract
	@param tin_dict <dict>: Dictionary to populate where [samplebasename][transcriptid] = tin_value
	@param key_index <int>: Index of the field to join multiple files
	@param parse_index <int>: Index of field of interest (i.e. TIN value)
	"""

	with open(file, 'r') as fh:
		header = next(fh).strip().split('\t')
		colid = header[key_index]
		file = os.path.basename(file) # Remove PATH
		sample = file.split(".star_rg_added.sorted.dmark.tin.xls")[0]

		for line in fh:
			linelist = line.strip().split('\t')
			tid = linelist[key_index]
			tinvalue = linelist[parse_index]
			if sample not in tin_dict:
				tin_dict[sample] = {}

			tin_dict[sample][tid] = tinvalue

	return colid, tin_dict




if __name__ == '__main__':

	# Get filenames to parse
	args = sys.argv
	files = sys.argv[1:]

	# Check if at least two files were provided
	if not len(args) >= 2:
		print("FATAL: Failed to provide more than one input file!")
		sys.exit("Usage:\n python {} *.tin.xls > combinedTIN.tsv".format(args[0]))

	# Populate tins with TINS values for all transcripts across all samples
	tins = {}
	for file in files:
		keycolname, tins = create(file, tins)

	df = pandas.DataFrame(tins)
	# Print dataframe to standard output
	df.to_csv(sys.stdout, sep="\t", header=True, index=True, index_label = keycolname)
