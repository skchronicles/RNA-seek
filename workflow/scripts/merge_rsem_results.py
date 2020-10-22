import os,sys
from os.path import join
from functools import reduce
import pandas as pd
if len(sys.argv) != 4:
	exit(os.path.basename(__file__)+":Incorrect number of arguments")
ens2genefile=sys.argv[1]
rsemgenesfolder=sys.argv[2]
rsemisoformsfolder=sys.argv[3]
annotations=pd.read_csv(ens2genefile,header=None,sep=" ",usecols=[0,2])
annotations.columns=["gene_id","GeneName"]
#collect gene results
files = list(filter(lambda x:"RSEM.genes.results" in x,os.listdir(rsemgenesfolder)))
print(files)
for cols in ["expected_count","TPM","FPKM"]:
    print(cols)
    dflist = list()
    dflist.append(annotations)
    for f in files:
        x=pd.read_csv(join(rsemgenesfolder,f),sep="\t",usecols=["gene_id",cols])
        samplename=f.split(".RSEM")[0]
        print(samplename)
        x.columns=["gene_id",samplename+"_"+cols]
        dflist.append(x)
    mergeddf=reduce(lambda a,b:pd.merge(a,b,how="outer",on="gene_id"),dflist)
    mergeddf.fillna('UNKNOWN',inplace=True)
    mergeddf=mergeddf.sort_values(by=['GeneName'])
    outfile=join(rsemgenesfolder,"RSEM.genes."+cols+".all_samples.txt")
    mergeddf.to_csv(outfile,sep="\t",index=False)
#collect isoform results
files = list(filter(lambda x:"RSEM.isoforms.results" in x,os.listdir(rsemisoformsfolder)))
print(files)
for cols in ["expected_count","TPM","FPKM"]:
    print(cols)
    dflist = list()
    for f in files:
        x=pd.read_csv(join(rsemisoformsfolder,f),sep="\t",usecols=["transcript_id","gene_id",cols])
        samplename=f.split(".RSEM")[0]
        print(samplename)
        x.columns=["transcript_id","gene_id",samplename+"_"+cols]
        dflist.append(x)
    mergeddf=reduce(lambda a,b:pd.merge(a,b,how="outer",on=["transcript_id","gene_id"]),dflist)
    print("reduced")
    mergeddf=pd.merge(annotations,mergeddf,on="gene_id")
    print("merged")
    mergeddf.fillna('UNKNOWN',inplace=True)
    mergeddf=mergeddf.sort_values(by=['GeneName'])
    print("sorted")
    outfile=join(rsemisoformsfolder,"RSEM.isoforms."+cols+".all_samples.txt")
    mergeddf.to_csv(outfile,sep="\t",index=False)
