from __future__ import print_function
import sys

# USAGE
# sys.argv[1] = genes.gtf

# Example
# $ python gtf2protein_coding_genes.py genes.gtf > protein_coding_genes.lst


def get_value(mykey, lookup):
	try:
		myvalue = lookup[mykey]
	except KeyError:
		myvalue = ''
	return myvalue.strip('"').strip("'")


def seperated(pairslist):
	for kv in pairslist:
		k = kv.split(' ')[0]
		v = " ".join(kv.split(' ')[1:]).rstrip(';')
		yield k,v


def get_id_and_type(last_column):
	pairs = {}
	kv_pairs_list = last_column.strip().split('; ')

	for k,v in seperated(kv_pairs_list):
		pairs[k] = v

	gene_id = get_value('gene_id', pairs)
	gene_type = get_value('gene_type', pairs)
	if not gene_type:
		# gene_type does not exist
		# default to using gene_biotype
		gene_type = get_value('gene_biotype', pairs)

	return gene_id, gene_type

if __name__ == '__main__':


	if len(sys.argv) != 2:
		print('Usage: python {} genes.gtf > protein_coding_genes.lst'.format(sys.argv[0]))
		print('\nError: failed to provide all positional arguments!', file=sys.stderr)
		sys.exit(1)

	protein_coding_genes = []
	with open(sys.argv[1]) as file:
		for line in file:
			if line.startswith('#'):
				# Skip over comments in header section
				continue

			linelist = line.strip().split("\t")
			if linelist[2]=="gene":
				# Get gene_id and gene_type
				gene_id, gene_type = get_id_and_type(last_column = linelist[-1])
				if gene_type=="protein_coding":
					protein_coding_genes.append(gene_id)
					print(gene_id,)
