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
		
genelist=[]	
chrs=[]
f=open("karyobed.bed",'w')
for i in open(sys.argv[1]).readlines():
	if i.startswith("#"):
		continue
	j=i.strip().split("\t")
	if j[2]=="gene":
		start=j[3]
		end=j[4]
		gene_name=get_gene_name(j[-1])
		if not gene_name in genelist:
			genelist.append(gene_name)
			outtxt=[]
			outtxt.append(j[0])
			outtxt.append(start)
			outtxt.append(end)
			outtxt.append(j[6])
			outtxt.append(gene_name)
			f.write("\t".join(outtxt)+"\n")
			chrs.append(j[0])
f.close()
chrs=list(set(chrs))
for c in chrs:
	f=open("karyobed."+c+".bed",'w')
	for i in open("karyobed.bed").readlines():
		j=i.strip().split("\t")
		if j[0]==c:
			f.write(i)
	f.close()
	
