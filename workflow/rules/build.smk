configfile:"config/build.yml"
GENOME=config["GENOME"].strip().replace(' ', '')
READLENGTHS=config["READLENGTHS"]
REFFA=config["REFFA"]
GTFFILE=config["GTFFILE"]
GTFVER=config["GTFVER"].strip().replace(' ', '')
OUTDIR=config["OUTDIR"]
SCRIPTSDIR=config["SCRIPTSDIR"]
workdir:OUTDIR


rule all:
	input:
		expand("rsemref/{genome}.transcripts.ump",genome=GENOME),
		"annotate.isoforms.txt",
		"annotate.genes.txt",
		"refFlat.txt",
		"geneinfo.bed",
		expand("STAR/2.7.0f/genes-{readlength}/SA",readlength=READLENGTHS),
		expand("{genome}.rRNA_interval_list",genome=GENOME),
		"karyoplot_gene_coordinates.txt",
		"qualimap_info.txt",
		"karyobeds/karyobed.bed",
		expand("{genome}_{gtfver}.json",genome=GENOME, gtfver=GTFVER)


rule rsem:
	input:
		fa=REFFA,
		gtf=GTFFILE
	params:
		genome=GENOME
	threads: 32
	output:
		"rsemref/{sample}.transcripts.ump"
	shell:'''
if [ ! -d rsemref ]; then mkdir rsemref; fi
cd rsemref
if [ ! -f {input.gtf} ]; then ln -s ../{input.gtf} .; fi
if [ ! -f {input.fa} ]; then ln -s ../{input.fa} .; fi
module load rsem/1.3.0
rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.fa} {params.genome}
rsem-generate-ngvector {params.genome}.transcripts.fa {params.genome}.transcripts
'''


rule annotate:
	input:
		gtf=GTFFILE
	output:
		"annotate.isoforms.txt",
		"annotate.genes.txt",
		"refFlat.txt",
		"genes.ref.bed",
		"geneinfo.bed"
	params:
		sdir=SCRIPTSDIR
	shell:'''
module load python/3.6
python {params.sdir}/get_gene_annotate.py {input.gtf} > annotate.genes.txt
python {params.sdir}/get_isoform_annotate.py {input.gtf} > annotate.isoforms.txt
module load ucsc/384
gtfToGenePred -ignoreGroupsWithoutExons {input.gtf} genes.genepred
genePredToBed genes.genepred genes.bed12
sort -k1,1 -k2,2n genes.bed12 > genes.ref.bed
python {params.sdir}/make_refFlat.py > refFlat.txt
python {params.sdir}/make_geneinfo.py > geneinfo.bed
'''


rule star:
	input:
		fa=REFFA,
		gtf=GTFFILE
	threads: 32
	output:
		SA="STAR/2.7.0f/genes-{readlength}/SA"
	shell:'''
# Creat Index for each readlength 
rl={wildcards.readlength}
rl=$((rl-1))

mkdir -p STAR/2.7.0f
cd STAR/2.7.0f
module load STAR/2.7.0f

STAR \
--runThreadN {threads} \
--runMode genomeGenerate \
--genomeDir ./genes-{wildcards.readlength} \
--genomeFastaFiles {input.fa} \
--sjdbGTFfile {input.gtf} \
--sjdbOverhang $rl \
--outTmpDir tmp_{wildcards.readlength} 
'''


rule rRNA_list:
	input:
		fa=REFFA,
		gtf=GTFFILE,
	output:
		expand("{genome}.rRNA_interval_list",genome=GENOME)
	params:
		genome=GENOME,
		sdir=SCRIPTSDIR
	shell:'''
module load samtools/1.9
module load python/3.6
python {params.sdir}/create_rRNA_intervals.py {input.fa} {input.gtf} {params.genome} > {params.genome}.rRNA_interval_list
'''


rule karyo_coord:
	input:
		gtf=GTFFILE
	output:
		"karyoplot_gene_coordinates.txt"
	params:
		sdir=SCRIPTSDIR
	shell:'''
module load python/3.6
python {params.sdir}/get_karyoplot_gene_coordinates.py {input.gtf} > karyoplot_gene_coordinates.txt
'''


rule qualimapinfo:
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		"qualimap_info.txt"
	params:
		sdir=SCRIPTSDIR
	shell:'''
module load python/2.7
python {params.sdir}/generate_qualimap_ref.py -g {input.gtf} -f {input.fa} -o {output} --ignore-strange-chrom 2> qualimap_error.log
'''


rule karyo_beds:
	input:
		gtf=GTFFILE
	output:
		"karyobeds/karyobed.bed"
	params:
		sdir=SCRIPTSDIR
	shell:'''
module load python/3.6
mkdir -p karyobeds
cd karyobeds
python {params.sdir}/get_karyoplot_beds.py {input.gtf}
'''


rule jsonmaker:
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		json="{sample}.json"
	params:
		workdir=OUTDIR,
		genome=GENOME
	run:
		import json
		outdir=params.workdir
		if not outdir.endswith("/"):
			outdir+="/"
		bigdict=dict()
		bigdict["references"]=dict()
		for i in ["exomeseq", "genomeseq", "rnaseq", "rnaseqvargerm", "ChIPseq"]:
			bigdict["references"][i]=dict()
		bigdict["references"]["rnaseq"]["GENOMEFILE"]=input.fa
		bigdict["references"]["rnaseq"]["GENOME"]=input.fa
		bigdict["references"]["rnaseq"]["GTFFILE"]=input.gtf
		bigdict["references"]["rnaseq"]["STARDIR"]=outdir+"STAR/2.7.0f/genes-"
		bigdict["references"]["rnaseq"]["STARREF"]=outdir+"STAR/2.7.0f/genes-"
		bigdict["references"]["rnaseq"]["ANNOTATE"]=outdir+"annotate.genes.txt"
		bigdict["references"]["rnaseq"]["ANNOTATEISOFORMS"]=outdir+"annotate.isoforms.txt"
		bigdict["references"]["rnaseq"]["REFFLAT"]=outdir+"refFlat.txt"
		bigdict["references"]["rnaseq"]["BEDREF"]=outdir+"genes.ref.bed"
		bigdict["references"]["rnaseq"]["GENEINFO"]=outdir+"geneinfo.bed"
		bigdict["references"]["rnaseq"]["QUALIMAP_INFO"]=outdir+"qualimap_info.txt"
		bigdict["references"]["rnaseq"]["KARYOBEDS"]=outdir+"karyobeds/"
		bigdict["references"]["rnaseq"]["KARYOPLOTER"]=outdir+"karyoplot_gene_coordinates.txt"
		bigdict["references"]["rnaseq"]["RSEMREF"]=outdir+"rsemref/"+params.genome
		bigdict["references"]["rnaseq"]["RRNALIST"]=outdir+params.genome+".rRNA_interval_list"
		bigdict["references"]["rnaseq"]["FASTQ_SCREEN_CONFIG"]="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf"
		bigdict["references"]["rnaseq"]["FASTAWITHADAPTERSETC"]="/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters_new.fa"
		bigdict["references"]["rnaseq"]["adapter.file"]="/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters.ngsqc.dat"
		bigdict["references"]["rnaseq"]["trimmomatic.adapters"]="/data/CCBR_Pipeliner/db/PipeDB/dev/adapters2.fa"
		bigdict["references"]["rnaseq"]["fastqc.adapters"]="/data/CCBR_Pipeliner/db/PipeDB/dev/fastqc.adapters"
		bigdict["references"]["rnaseq"]["ORGANISM"]="CUSTOM"
		with open(output.json, 'w') as fp:
			json.dump(bigdict, fp, indent=4)
