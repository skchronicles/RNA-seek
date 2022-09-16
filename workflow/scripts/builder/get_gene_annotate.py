import sys;
unknown="\"Unknown\";"

for i in list(filter(lambda x:x[2]=="gene",filter(lambda x:not x[0].startswith("#"),list(map(lambda x:x.strip().split("\t"),open(sys.argv[1]).readlines()))))):
#for i in list(map(lambda x:x.strip().split("\t")[8],open(sys.argv[1]).readlines())):
	gene_id=""
	j=i[8].split()
	gene_id=unknown
	gene_name=unknown
	gene_biotype=unknown
	for k in list(range(0,len(j)-1,2)):
		if j[k]=="gene_id":
			gene_id=j[k+1]
		elif j[k]=="gene_name":
			gene_name=j[k+1]
		elif j[k]=="gene_biotype":
			gene_biotype=j[k+1]
		elif j[k]=="gene_type":
			gene_biotype=j[k+1]
	if gene_name == unknown and gene_id != unknown :
		gene_name=gene_id
	s="%s  %s  %s"%(gene_id[:-1],gene_name[:-1],gene_biotype[:-1])
	print(s)
	
