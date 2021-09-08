#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from __future__ import print_function
import pandas as pd
import xlsxwriter
import sys
import os


def reader(filename, subset=[], skip='#', **kwargs):
	"""Reads in an MAF-like file as a dataframe. Determines the 
	correct handler for reading in a given MAF file. Supports reading
	in TSV files (.tsv, .txt, .text, .vcf, or .maf), CSV files (.csv), 
	and excel files (.xls, .xlsx, .xlsm, .xlsb, .odf, .ods, .odt ). 
	The subset option allows a users to only select a few columns 
	given a list of column names.
	@param filename <str>:
		Path of an MAF-like file to read and parse
	@param subset list[<str>]:
		List of column names which can be used to subset the df
	@param skip <str>:
		Skips over line starting with this character
	@params kwargs <read_excel()>
		Key words to modify pandas.read_excel() function behavior
	@return <pandas dataframe>:
		dataframe with spreadsheet contents
	"""
	# Get file extension
	extension = os.path.splitext(filename)[-1].lower()

	# Assign a handler to read in the file
	if extension in ['.xls', '.xlsx', '.xlsm', '.xlsb', '.odf', '.ods', '.odt']:
		# Read in as an excel file
		return excel(filename, subset, skip, **kwargs)
	elif extension in ['.csv']:
		# Read in as an CSV file
		return csv(filename, subset, skip, **kwargs)
	else:
		# Default to reading in as an TSV file
		# Tab is the normal delimeter for MAF or VCF files
		# MAF files usually have one of the following
		# extensions: '.tsv', '.txt', '.text', '.vcf', '.maf'
		return tsv(filename, subset, skip, **kwargs)


def excel(filename, subset=[], skip='#', **kwargs):
	"""Reads in an excel file as a dataframe. The subset option
	allows a users to only select a few columns given a list of 
	column names.
	@param filename <str>:
		Path of an EXCEL file to read and parse
	@param subset list[<str>]:
		List of column names which can be used to subset the df
	@param skip <str>:
		Skips over line starting with this character
	@params kwargs <read_excel()>
		Key words to modify pandas.read_excel() function behavior
	@return <pandas dataframe>:
		dataframe with spreadsheet contents
	"""
	if subset:
		return pd.read_excel(filename, comment=skip, **kwargs)[subset]

	return pd.read_excel(filename, comment=skip, **kwargs)


def tsv(filename, subset=[], skip='#', **kwargs):
	"""Reads in an TSV file as a dataframe. The subset option
	allows a users to only select a few columns given a list of 
	column names.
	@param filename <str>:
		Path of an TSV file to read and parse
	@param subset list[<str>]:
		List of column names which can be used to subset the df
	@param skip <str>:
		Skips over line starting with this character
	@params kwargs <read_excel()>
		Key words to modify pandas.read_excel() function behavior
	@return <pandas dataframe>:
		dataframe with spreadsheet contents
	"""
	if subset:
		return pd.read_table(filename, comment=skip, **kwargs)[subset]

	return pd.read_table(filename, comment=skip, **kwargs)


def csv(filename, subset=[], skip='#', **kwargs):
	"""Reads in an CSV file as a dataframe. The subset option
	allows a users to only select a few columns given a list of 
	column names.
	@param filename <str>:
 		Path of an CSV file to read and parse
	@param subset list[<str>]:
		List of column names which can be used to subset the df
	@param skip <str>:
		Skips over line starting with this character
	@params kwargs <read_excel()>
		Key words to modify pandas.read_excel() function behavior
	@return <pandas dataframe>:
		dataframe with spreadsheet contents
	"""
	if subset:
		return pd.read_csv(filename, comment=skip, **kwargs)[subset]

	return pd.read_csv(filename, comment=skip, **kwargs)


def excel_writer(files, spreadsheet= 'test.xlsx'):
	"""Takes a list of files and creates one excel spreadsheet.
	Each file will becomes a sheet in the spreadsheet where the 
	name of the sheet is the basename of the file with the extension
	removed.
	@param files list[<str>]:
    	List of files to merge into one execl file
	@param spreadsheet <str>:
    	Output filename of the spreadsheet
	"""

	writer = pd.ExcelWriter(spreadsheet, engine='xlsxwriter')
	
	# Create a spreadsheet from the contents of each file
	for file in files:
		print('Reading in {}'.format(file))
		df = reader(file)
		sheet = os.path.splitext(os.path.basename(file))[0]
		try:
			# Sheet name cannot exceed 31 characters in length
			df.to_excel(writer, sheet_name = sheet, index = False, freeze_panes = (1,0))
		except xlsxwriter.exceptions.InvalidWorksheetName as e:
			df.to_excel(writer, sheet_name = sheet[:31], index = False,freeze_panes = (1,0))

	writer.save()


if __name__ == '__main__':

	# List of file to convert into an excel file
	files = sys.argv[1:-1]
	# Output file name 
	outfh = sys.argv[-1]

	excel_writer(files, outfh)
