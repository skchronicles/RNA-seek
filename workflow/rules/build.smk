from os.path import join, basename

# Global Workflow variables
configfile:join("config","build.yml")
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
		# expand("STAR/2.7.6a/genes-{readlength}/SA",readlength=READLENGTHS),
		expand("{genome}_{gtfver}.json",genome=GENOME, gtfver=GTFVER),
		expand("{genome}.rRNA_interval_list",genome=GENOME),
		expand("rsemref/{genome}.transcripts.ump",genome=GENOME),
		"STAR/2.7.6a/genome/SA",
		"annotate.isoforms.txt",
		"annotate.genes.txt",
		"refFlat.txt",
		"geneinfo.bed",
		"karyoplot_gene_coordinates.txt",
		"qualimap_info.txt",
		"karyobeds/karyobed.bed",
		"transcripts.protein_coding_only.bed12"


rule rsem:
    """
    Builds reference files for rsem to quantify gene and isoform counts.
	@Input:
		Genomic FASTA file
		Annotation file in GTF format
	@Output:
        Generates 'reference_name.grp', 'reference_name.ti',
    	'reference_name.transcripts.fa', 'reference_name.seq',
    	'reference_name.chrlist', 'reference_name.idx.fa',
    	'reference_name.n2g.idx.fa', 'reference_name.grp',
		'reference_name.ti', 'reference_name.seq', and
		'reference_name.chrlist' which are used by RSEM internally.
		'reference_name.transcripts.fa' contains the extracted reference
		transcripts in Multi-FASTA format.
	"""
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		join("rsemref","{genome}.transcripts.ump")
	threads: 32
	params:
		rname='bl:rsem',
		genome=GENOME,
		prefix=join("rsemref", GENOME)
	container: "docker://nciccbr/ccbr_rsem_1.3.3:v1.0"
	shell: """
	rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.fa} {params.prefix}
	rsem-generate-ngvector {params.prefix}.transcripts.fa {params.prefix}.transcripts
	"""


rule annotate:
	"""
	Builds reference files for rsem_merge (annotate.genes.txt),
	stats (refFlat.txt), rseqc inner_distance (genes.ref.bed),
	rseqc read_distribution (genes.ref.bed), rseqc infer_experiment
	(genes.ref.bed).
	@Input:
		Annotation file in GTF format
	@Output:
		Generates 'annotate.isoforms.txt', 'annotate.genes.txt',
		'refFlat.txt', 'genes.ref.bed', 'geneinfo.bed'
	"""
	input:
		gtf=GTFFILE
	output:
		"annotate.isoforms.txt",
		"annotate.genes.txt",
		"refFlat.txt",
		"genes.ref.bed",
		"geneinfo.bed"
	params:
		rname='bl:annotate',
		get_gene=join(SCRIPTSDIR, "get_gene_annotate.py"),
		get_isoform=join(SCRIPTSDIR, "get_isoform_annotate.py"),
		make_refFlat=join(SCRIPTSDIR, "make_refFlat.py"),
		make_geneinfo=join(SCRIPTSDIR, "make_geneinfo.py")
	container: "docker://nciccbr/ccbr_build_rnaseq:v0.0.1"
	shell: """
	python3 {params.get_gene} {input.gtf} > annotate.genes.txt
	python3 {params.get_isoform} {input.gtf} > annotate.isoforms.txt
	gtfToGenePred -ignoreGroupsWithoutExons {input.gtf} genes.genepred
	genePredToBed genes.genepred genes.bed12
	sort -k1,1 -k2,2n genes.bed12 > genes.ref.bed
	python3 {params.make_refFlat} > refFlat.txt
	python3 {params.make_geneinfo} {input.gtf} > geneinfo.bed
	"""


rule star_rl:
    """
    Builds STAR Index to align reads against reference genome for a defined
	read length. Optimal readlength for sjdbOverhang = max(ReadLength) - 1.
	For most applications, a readlength of 100 works well.
    @Input:
    	Genomic FASTA file
		Annotation file in GTF format
    @Output:
        Generates 'chrName.txt', 'chrLength.txt', 'chrStart.txt',
		'chrNameLength.txt', 'exonGeTrInfo.tab', 'geneInfo.tab',
		'transcriptInfo.tab', 'exonInfo.tab', 'sjdbList.fromGTF.out.tab',
		'sjdbInfo.txt', 'sjdbList.out.tab', 'genomeParameters.txt',
		'Genome', 'SA', 'SAindex' which are used by STAR internally
		for a list of provided readlengths.
	"""
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		join("STAR", "2.7.6a", "genes-{readlength}", "SA")
	threads: 32
	params:
		rname='bl:star_rl',
	container: "docker://nciccbr/ccbr_arriba_2.0.0:v0.0.1"
	shell: """
	# Create Index for read length
	rl=$(({wildcards.readlength}-1))

	STAR \
		--runThreadN {threads} \
		--runMode genomeGenerate \
		--genomeDir STAR/2.7.6a/genes-{wildcards.readlength} \
		--genomeFastaFiles {input.fa} \
		--sjdbGTFfile {input.gtf} \
		--sjdbOverhang $rl \
		--outFileNamePrefix STAR/2.7.6a/build_{wildcards.readlength}_ \
		--outTmpDir /lscratch/$SLURM_JOB_ID/tmp_{wildcards.readlength}
	"""


rule star_genome:
    """
    Builds STAR Index to align reads against reference genome without the
	GTF or readlength provided. This index only contain information pertaining
	to the asssembly or reference genome in FASTA format. This indice represents a
	base index from which processing annotations from a GTF and insert junctions
	on the fly. This has the advantage of saving diskspace as an index will not be
	created for a list of predefined readlengths. This rule replaces star_rl above.
    @Input:
    	Genomic FASTA file
    @Output:
        Generates 'chrName.txt', 'chrLength.txt', 'chrStart.txt',
		'chrNameLength.txt', 'exonGeTrInfo.tab', 'geneInfo.tab',
		'transcriptInfo.tab', 'exonInfo.tab', 'sjdbList.fromGTF.out.tab',
		'sjdbInfo.txt', 'sjdbList.out.tab', 'genomeParameters.txt',
		'Genome', 'SA', 'SAindex' which are used by STAR internally.
	"""
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		join("STAR", "2.7.6a", "genome", "SA")
	threads: 32
	params:
		rname='bl:star_genome',
	container: "docker://nciccbr/ccbr_arriba_2.0.0:v0.0.1"
	shell: """
	STAR \
		--runThreadN {threads} \
		--runMode genomeGenerate \
		--genomeDir STAR/2.7.6a/genome \
		--genomeFastaFiles {input.fa} \
		--outFileNamePrefix STAR/2.7.6a/build_genome_ \
		--outTmpDir /lscratch/$SLURM_JOB_ID/tmp_genome
	"""


rule rRNA_list:
    """
    Builds ribosomal RNA list reference file for picard CollectRnaMetrics to
	estimate rRNA content or abundance. rRNA can make up a significant proportion
	of an RNA-seq library if not properly depleted either through poly-selection
	or ribosomal depletion. Samples with very high rRNA content could signal an issue
	occured with library prepartion.
    @Input:
    	Genomic FASTA file
		Annotation file in GTF format
    @Output:
        Generates '{genome}.rRNA_interval_list' which is used by picard internally.
	"""
	input:
		fa=REFFA,
		gtf=GTFFILE,
	output:
		"{genome}.rRNA_interval_list"
	params:
		rname='bl:rRNA_list',
		genome=GENOME,
		create_rRNA=join(SCRIPTSDIR, "create_rRNA_intervals.py")
	container: "docker://nciccbr/ccbr_build_rnaseq:v0.0.1"
	shell: """
	python3 {params.create_rRNA} \
		{input.fa} \
		{input.gtf} \
		{params.genome} > {params.genome}.rRNA_interval_list
	"""


rule karyo_coord:
    """
    Builds karyoplot reference file used by limma, edgeR, DESeq2 DEG reports.
    @Input:
		Annotation file in GTF format
    @Output:
        Generates 'karyoplot_gene_coordinates.txt' which is used by DE reports internally.
	"""
	input:
		gtf=GTFFILE
	output:
		"karyoplot_gene_coordinates.txt"
	params:
		rname='bl:karyo_coord',
		get_karyoplot=join(SCRIPTSDIR, "get_karyoplot_gene_coordinates.py")
	container: "docker://nciccbr/ccbr_build_rnaseq:v0.0.1"
	shell: """
	python3 {params.get_karyoplot} {input.gtf} > karyoplot_gene_coordinates.txt
	"""


rule karyo_beds:
    """
    Builds karyoplot reference file used by limma, edgeR, DESeq2 DEG reports.
    @Input:
		Annotation file in GTF format
    @Output:
        Generates 'karyobeds/karyobed.bed' which is used by DE reports internally.
	"""
	input:
		gtf=GTFFILE
	output:
		join("karyobeds", "karyobed.bed")
	params:
		rname='bl:karyo_bed',
		get_karyoplot=join(SCRIPTSDIR, "get_karyoplot_beds.py")
	container: "docker://nciccbr/ccbr_build_rnaseq:v0.0.1"
	shell: """
	mkdir -p karyobeds && cd karyobeds
	python3 {params.get_karyoplot} {input.gtf}
	"""


rule tin_ref:
    """
    Builds RSeQC tin.py reference file containing all canocical protein coding genes.
    @Input:
		Annotation file in GTF format
    @Output:
        Generates 'transcripts.protein_coding_only.bed12' which is used by RSeQC tin.py internally.
	"""
	input:
		gtf=GTFFILE
	output:
		"transcripts.protein_coding_only.bed12"
	params:
		rname='bl:tin_ref',
		gtf2protein=join(SCRIPTSDIR, "gtf2protein_coding_genes.py"),
		gene2transcripts=join(SCRIPTSDIR, "gene2transcripts_add_length.py"),
	container: "docker://nciccbr/ccbr_build_rnaseq:v0.0.1"
	shell: """
	python {params.gtf2protein} {input.gtf} > protein_coding_genes.lst
	gtfToGenePred -genePredExt {input.gtf} genes.gtf.genePred
	genePredToBed genes.gtf.genePred genes.gtf.genePred.bed

	awk -F '\\t' -v OFS='\\t' '{{print $12,$1}}' genes.gtf.genePred \
		| sort -k1,1n > gene2transcripts

	while read gene; do
		grep "${{gene}}" gene2transcripts;
	done < protein_coding_genes.lst > gene2transcripts.protein_coding_only

	python {params.gene2transcripts} gene2transcripts.protein_coding_only genes.gtf.genePred.bed \
		> gene2transcripts.protein_coding_only.with_len

	sort -k1,1 -k3,3nr gene2transcripts.protein_coding_only.with_len | \
		awk -F '\\t' '{{if (!seen[$1]) {{seen[$1]++; print $2}}}}' > protein_coding_only.txt

	while read transcript; do
		grep -m1 "${{transcript}}" genes.gtf.genePred.bed;
	done < <(awk -F '.' '{{print $1}}' protein_coding_only.txt) > transcripts.protein_coding_only.bed12

	rm -f "protein_coding_genes.lst" "genes.gtf.genePred" "genes.gtf.genePred.bed" \
		"gene2transcripts" "gene2transcripts.protein_coding_only" "protein_coding_only.txt" \
		"gene2transcripts.protein_coding_only.with_len" "protein_coding_only.txt"
	"""


rule qualimapinfo:
    """
    Builds QualiMap2 reference file used by Counts QC. It is a four-column
	tab-delimited text file, with the features names or IDs in the first column,
	the group (e.g. the biotype from Ensembl database) in the second column,
	feature length in the third column and feature GC-content in the last column.
    @Input:
    	Genomic FASTA file
		Annotation file in GTF format
    @Output:
        Generates 'qualimap_info.txt' which is used by QualiMap2 internally.
	"""
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		"qualimap_info.txt"
	params:
		rname='bl:qualimapinfo',
		generate_qualimap=join(SCRIPTSDIR, "generate_qualimap_ref.py")
	container: "docker://nciccbr/ccbr_build_rnaseq:v0.0.1"
	shell: """
	python3 {params.generate_qualimap} \
		-g {input.gtf} \
		-f {input.fa} \
		-o {output} \
		--ignore-strange-chrom 2> qualimap_error.log
	"""


rule jsonmaker:
    """
    Builds reference genome reference JSON file. This is a config file for the
	RNA-seek pipeline that defines each location to a given reference file.
    @Input:
	    Genomic FASTA file
		Annotation file in GTF format
    @Output:
        Generates '{genome}.json' which is used RNA-seek pipeline internally.
	"""
	input:
		fa=REFFA,
		gtf=GTFFILE
	output:
		json="{genome}.json"
	params:
		rname='bl:jsonmaker',
		workdir=OUTDIR,
		genome=GENOME
	run:
		import json
		outdir=params.workdir
		if not outdir.endswith("/"):
			outdir+="/"
		refdict = {}
		refdict["references"] = {}
		refdict["references"]["rnaseq"] = {}
		refdict["references"]["rnaseq"]["GENOMEFILE"] = input.fa
		refdict["references"]["rnaseq"]["GENOME"] = input.fa
		refdict["references"]["rnaseq"]["GTFFILE"] = input.gtf
		refdict["references"]["rnaseq"]["GENOME_STARDIR"] = outdir+"STAR/2.7.6a/genome"
		refdict["references"]["rnaseq"]["ANNOTATE"] = outdir+"annotate.genes.txt"
		refdict["references"]["rnaseq"]["ANNOTATEISOFORMS"] = outdir+"annotate.isoforms.txt"
		refdict["references"]["rnaseq"]["REFFLAT"] = outdir+"refFlat.txt"
		refdict["references"]["rnaseq"]["BEDREF"] = outdir+"genes.ref.bed"
		refdict["references"]["rnaseq"]["GENEINFO"] = outdir+"geneinfo.bed"
		refdict["references"]["rnaseq"]["QUALIMAP_INFO"] = outdir+"qualimap_info.txt"
		refdict["references"]["rnaseq"]["KARYOBEDS"] = outdir+"karyobeds/"
		refdict["references"]["rnaseq"]["KARYOPLOTER"] = outdir+"karyoplot_gene_coordinates.txt"
		refdict["references"]["rnaseq"]["RSEMREF"] = outdir+"rsemref/"+params.genome
		refdict["references"]["rnaseq"]["RRNALIST"] = outdir+params.genome+".rRNA_interval_list"
		refdict["references"]["rnaseq"]["ORGANISM"] = wildcards.genome
		refdict["references"]["rnaseq"]["TINREF"] = outdir+"transcripts.protein_coding_only.bed12"

		# Try to infer which Arriba reference files to add a user defined reference genome
		if 'hg19' in params.genome.lower() or \
		'hs37d' in params.genome.lower() or \
		'grch37' in params.genome.lower():
			refdict["references"]["rnaseq"]["FUSIONBLACKLIST"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/blacklist_hg19_hs37d5_GRCh37_v2.0.0.tsv.gz"
			refdict["references"]["rnaseq"]["FUSIONCYTOBAND"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/cytobands_hg19_hs37d5_GRCh37_v2.0.0.tsv"
			refdict["references"]["rnaseq"]["FUSIONPROTDOMAIN"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.0.0.gff3"
		elif 'hg38' in params.genome.lower() or \
		'hs38d' in params.genome.lower() or \
		'grch38' in params.genome.lower():
			refdict["references"]["rnaseq"]["FUSIONBLACKLIST"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/blacklist_hg38_GRCh38_v2.0.0.tsv.gz"
			refdict["references"]["rnaseq"]["FUSIONCYTOBAND"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/cytobands_hg38_GRCh38_v2.0.0.tsv"
			refdict["references"]["rnaseq"]["FUSIONPROTDOMAIN"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/protein_domains_hg38_GRCh38_v2.0.0.gff3"
		elif 'mm10' in params.genome.lower() or \
		'grcm38' in params.genome.lower():
			refdict["references"]["rnaseq"]["FUSIONBLACKLIST"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/blacklist_mm10_GRCm38_v2.0.0.tsv.gz"
			refdict["references"]["rnaseq"]["FUSIONCYTOBAND"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/cytobands_mm10_GRCm38_v2.0.0.tsv"
			refdict["references"]["rnaseq"]["FUSIONPROTDOMAIN"] = \
			"s3://nciccbr/Resources/RNA-seq/arriba/protein_domains_mm10_GRCm38_v2.0.0.gff3"

		with open(output.json, 'w') as fp:
			json.dump(refdict, fp, indent=4)
