annotate_genes={a[0]:a for a in map(lambda x:x.strip().replace('"','').split("  "),open("annotate.genes.txt").readlines())}

for l in list(filter(lambda x:x[2]=="gene",filter(lambda x:not x[0].startswith("#"),list(map(lambda x:x.strip().split("\t"),open("genes.gtf").readlines()))))):
	newl=[]
	newl.append(l[0])
	newl.append(l[3])
	newl.append(l[4])
	newl.append(l[6])
	col9=l[8].split(" ")
	gene_id_index=col9.index("gene_id")
	gene_id=col9[gene_id_index+1].strip(";").strip("\"")
	newl.append(gene_id)
	try:
		#newl.append(annotate_genes["\""+gene_id+"\""][2].strip("\""))
		newl.append(annotate_genes[gene_id][2].strip("\""))
	except IndexError:
		print(gene_id)
		exit()
	newl.append(annotate_genes[gene_id][1].strip("\""))
	print("\t".join(newl))
