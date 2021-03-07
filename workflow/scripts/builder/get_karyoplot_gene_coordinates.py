import sys
def get_gene_name(j):
	searchfor="gene_name"
	if not searchfor in j:
		searchfor="gene_id"
	k=j.split()
	ind=-1
	for i,l in enumerate(k):
		if l==searchfor:
			ind=i+1
			break
	m=k[ind].split("\"")[1]
	return m
		
print("chr","coord","gene_name","strand",sep="\t")
genelist=[]	
for i in open(sys.argv[1]).readlines():
	if i.startswith("#"):
		continue
	j=i.strip().split("\t")
	if j[2]=="gene":
		coord=int((int(j[3])+int(j[4]))*0.5)
		gene_name=get_gene_name(j[-1])
		if not gene_name in genelist:
			genelist.append(gene_name)
			print(j[0],coord,gene_name,j[6],sep="\t")
