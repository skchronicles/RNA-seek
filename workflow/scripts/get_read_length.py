import sys,os,zipfile,glob

def get_max_read_length(f):
	returnvalues=[]
	for read in ["1","2"]:
		zipf=zipfile.ZipFile(f)
		foldername=os.path.splitext(os.path.basename(f))[0]
		fastqc_data=zipf.open(foldername+"/fastqc_data.txt")
		line=list(filter(lambda x:x.startswith(b'Sequence length'),list(map(lambda x:x.strip(),fastqc_data.readlines()))))[0]
		if b'-' in line:
			returnvalues.append(int(line.split(b'\t')[1].split(b'-')[1]))
		else:
			returnvalues.append(int(line.split(b'\t')[1]))
	return max(returnvalues)


basefolder=sys.argv[1] # QC or rawQC folder
fastqs=glob.glob(basefolder+"/*fastqc.zip")
read_lengths=[]
for f in fastqs:
    read_lengths.append(get_max_read_length(f))
print(max(read_lengths))
