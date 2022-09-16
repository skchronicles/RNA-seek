import sys;
unknown="\"Unknown\";"

for i in list(filter(lambda x:x[2]=="transcript",filter(lambda x:not x[0].startswith("#"),list(map(lambda x:x.strip().split("\t"),open(sys.argv[1]).readlines()))))):
#for i in list(map(lambda x:x.strip().split("\t")[8],open(sys.argv[1]).readlines())):
	gene_id=""
	j=i[8].split()
	transcript_id=unknown
	gene_id=unknown
	gene_name=unknown
	transcript_name=unknown
	for k in list(range(0,len(j)-1,2)):
		if j[k]=="transcript_id":
			transcript_id=j[k+1]
		elif j[k]=="gene_id":
			gene_id=j[k+1]
		elif j[k]=="transcript_name":
			transcript_name=j[k+1]
		elif j[k]=="gene_name":
			gene_name=j[k+1]
	if transcript_name == unknown and transcript_id != unknown :
		transcript_name=transcript_id
	if gene_name == unknown and gene_id != unknown :
		gene_name=gene_id
	s="%s  %s  %s  %s"%(transcript_id[:-1],transcript_name[:-1],gene_id[:-1],gene_name[:-1])
	print(s)
	
