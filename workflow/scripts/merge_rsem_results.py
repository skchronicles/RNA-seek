from __future__ import print_function
from os.path import join
from functools import reduce
import os,sys
import pandas as pd


def Counts(fpattern, searchpath, anno, ftype, mycols):
	"""
	Get each samples FPKM vaules from RSEMs *.RSEM.genes.results and *.RSEM.isoform.results
	"""
	# Collect RSEM Results
	files = list(filter(lambda x: fpattern in x,os.listdir(searchpath)))
	# print(files)
	for col in ["expected_count","TPM","FPKM"]:
	    dflist = []
	    for f in files:
	        x=pd.read_csv(join(searchpath,f),sep="\t",usecols=mycols+[col])
	        samplename=f.split(".RSEM")[0]
	        x.columns=mycols+[samplename]
	        dflist.append(x)

	    mergeddf=reduce(lambda a,b:pd.merge(a,b,how="outer",on=mycols),dflist)
	    mergeddf=pd.merge(anno,mergeddf,on="gene_id")
	    mergeddf.fillna('UNKNOWN',inplace=True)
	    mergeddf=mergeddf.sort_values(by=['GeneName'])
	    outfile=join(searchpath, "RSEM." + ftype + "." + col + ".all_samples.txt")
	    mergeddf.to_csv(outfile,sep="\t",index=False)


if __name__ == '__main__':
	# Parse Args
	if len(sys.argv) != 4:
		print(os.path.basename(__file__)+ ": Fail to provide all required arguments!")
		exit("USAGE: python {} /path/to/annotate.genes.txt /path/to/rsem/genecounts/ /path/to/isoformcounts/".format(sys.argv[0]))

	ens2genefile = sys.argv[1]
	rsemgenesfolder = sys.argv[2]
	rsemisoformsfolder = sys.argv[3]

	annotations=pd.read_csv(ens2genefile,header=None,sep=" ",usecols=[0,2])
	annotations.columns=["gene_id","GeneName"]

	Counts(fpattern = "RSEM.genes.results", searchpath = rsemgenesfolder, anno = annotations, ftype = "genes", mycols = ["gene_id"])
	Counts(fpattern = "RSEM.isoforms.results", searchpath = rsemisoformsfolder, anno = annotations, ftype = "isoforms", mycols = ["transcript_id","gene_id"])
