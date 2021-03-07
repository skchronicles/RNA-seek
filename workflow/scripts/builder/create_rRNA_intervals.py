import sys,os,pysam
fa=sys.argv[1]
genomename=sys.argv[3]
gtf=sys.argv[2]
if not os.path.exists(fa+".fai"):
	pysam.faidx(fa)
unknown="\"Unknown\";"
#out=open(genomename+".rRNA_interval_list",'w')
for f in open(fa+".fai").readlines():
	f=f.strip().split("\t")
#	out.write("@SQ\tSN:%s\tLN:%s\tAS:%s\n"%(f[0],f[1],genomename))
	print("@SQ\tSN:%s\tLN:%s\tAS:%s"%(f[0],f[1],genomename))
	

for i in list(filter(lambda x:x[2]=="gene",filter(lambda x:not x[0].startswith("#"),list(map(lambda x:x.strip().split("\t"),open(gtf).readlines()))))):
	gene_id=""
	j=i[8].split()
	gene_id=unknown
	gene_name=unknown
	gene_biotype=unknown
	for k in list(range(0,len(j)-1,2)):
		if j[k]=="gene_id":
			gene_id=j[k+1][1:-2]
		elif j[k]=="gene_name":
			gene_name=j[k+1][1:-2]
		elif j[k]=="gene_biotype":
			gene_biotype=j[k+1][1:-2]
		elif j[k]=="gene_type":
			gene_biotype=j[k+1][1:-2]
	if gene_biotype=="rRNA":
		#out.write("%s\t%s\t%s\t%s\t%s\n"%(i[0],i[3],i[4],i[6],gene_id))
		print("%s\t%s\t%s\t%s\t%s"%(i[0],i[3],i[4],i[6],gene_id))
#out.close()
