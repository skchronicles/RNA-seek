from __future__ import print_function
import sys

# USAGE:
# sys.argv[1] = gene2transcripts.protein_coding_only
# sys.argv[2] = genes.gtf.genePred.bed

# Example:
# $ python gene2transcripts_add_length.py gene2transcripts.protein_coding_only genes.gtf.genePred.bed > gene2transcripts.protein_coding_only.with_len

def get_len(s):
	lengths=[int(num) for num in s.strip().rstrip(',').split(",")]
	return sum(lengths)


if __name__ == '__main__':

	if len(sys.argv) != 3:
		print('Usage: python {} gene2transcripts.protein_coding_only genes.gtf.genePred.bed'.format(sys.argv[0]))
		print('\nError: failed to provide all positional arguments!', file=sys.stderr)
		sys.exit(1)

	transcript2length = {}

	# Read in genes.gtf.genePred.bed
	with open(sys.argv[2]) as in_file:
		for line in in_file:
			linelist = line.strip().split("\t")
			transcript2length[linelist[3]] = get_len(linelist[10])

	# Read in gene2transcripts.protein_coding_only
	with open(sys.argv[1]) as in_file:
		for line in in_file:
			linelist = line.strip().split("\t")
			linelist.append(str(transcript2length[linelist[1]]))
			print("\t".join(linelist))
