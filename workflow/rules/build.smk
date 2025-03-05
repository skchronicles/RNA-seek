from os.path import join, basename
import json


# Helper Functions
def allocated(resource, rule, lookup, default="__default__"):
    """Pulls resource information for a given rule. If a rule does not have any information 
    for a given resource type, then it will pull from the default. Information is pulled from
    definitions in the cluster.json (which is used a job submission). This ensures that any 
    resources used at runtime mirror the resources that were allocated.
    :param resource <str>: resource type to look in cluster.json (i.e. threads, mem, time, gres)
    :param rule <str>: rule to lookup its information
    :param lookup <dict>: Lookup containing allocation information (i.e. cluster.json)
    :param default <str>: default information to use if rule information cannot be found
    :return allocation <str>: 
        allocation information for a given resource type for a given rule
    """

    try: 
        # Try to get allocation information
        # for a given rule
        allocation = lookup[rule][resource]
    except KeyError:
        # Use default allocation information
        allocation = lookup[default][resource]
    
    return allocation


def provided(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    """   
    if not str_bool(condition):
        # If condition is False, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []
    return samplelist


def str_bool(s):
    """Converts a string to boolean. It is dangerous to try to
    typecast a string into a boolean value using the built-in 
    `bool()` function. This function avoids any issues that can
    arise when using `bool()`. 
    Example:
      boolean('True') returns True
      boolean('False') returns False
      boolean('asdas') raises TypeError
    """
    if not s:
        return False
    val = s.lower()
    if val in ['true', '1', 'y', 'yes']:
        return True
    elif val in ['false', '0', 'n', 'no', '', 'none']:
        return False
    else:
        # Provided a string or path
        return s


# Global Workflow variables
configfile:join("config","build.yml")
GENOME=config["GENOME"].strip().replace(' ', '')
READLENGTHS=config["READLENGTHS"]
REFFA=config["REFFA"]
GTFFILE=config["GTFFILE"]
GTFVER=config["GTFVER"].strip().replace(' ', '')
OUTDIR=config["OUTDIR"]
SCRIPTSDIR=config["SCRIPTSDIR"]
tmpdir=config["TMP_DIR"]
workdir:OUTDIR

# Read in resource information,
# containing information about 
# threads, mem, walltimes, etc.
# TODO: Add handler for when the
# mode is set to local.
with open(join(OUTDIR, 'resources', 'build_cluster.json')) as fh:
    cluster = json.load(fh)

# Ensures backwards compatibility 
try:
    SMALL_GENOME=config["SMALL_GENOME"]
except KeyError:
    SMALL_GENOME="False"

try:
    SHARED_PATH=config["SHARED_RESOURCES"]
except KeyError:
    SHARED_PATH="None"


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
        "transcripts.protein_coding_only.bed12",
        # FastQ Screen P1 Reference files,
        # conditional runs with --shared-resources option
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "hg19", "hg19.{ext}.bt2"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "mm9", "mm9.{ext}.bt2"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "Virus", "virus.{ext}.bt2"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "hg19", "hg19.rev.{ext}.bt2"), ext=["1", "2"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "mm9", "mm9.rev.{ext}.bt2"), ext=["1", "2"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "Virus", "virus.rev.{ext}.bt2"), ext=["1", "2"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "Bacteria", "bacteria.{ext}.bt2l"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "Fungi", "fungi.{ext}.bt2l"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "Bacteria", "bacteria.rev.{ext}.bt2l"), ext=["1", "2"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "Fungi", "fungi.rev.{ext}.bt2l"), ext=["1", "2"]), SHARED_PATH),
        # FastQ Screen P2 Reference files,
        # conditional runs with --shared-resources option
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "UniVec_vectors", "UniVec_vectors.{ext}.bt2"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "rRNA", "rRNA.{ext}.bt2"), ext=["1", "2", "3", "4"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "UniVec_vectors", "UniVec_vectors.rev.{ext}.bt2"), ext=["1", "2"]), SHARED_PATH),
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "rRNA", "rRNA.rev.{ext}.bt2"), ext=["1", "2"]), SHARED_PATH),
        # FastQ Screeen confs,
        # conditional runs with --shared-resources option
        provided(expand(join(SHARED_PATH, "fastq_screen_db", "fastq_screen_p{ext}.conf"), ext=["1", "2"]), SHARED_PATH),
        # Kraken2 Database,
        # conditional runs with --shared-resources option
        provided(expand(join(SHARED_PATH, "k2_pluspf_20241228", "{ref}.k2d"), ref=["hash", "opts", "taxo"]), SHARED_PATH),
        # Arriba reference files
        # conditional runs with --shared-resources option
        provided(expand(join(SHARED_PATH, "arriba", "blacklist_{ref}_v2.0.0.tsv.gz"), ref=["hg38_GRCh38", "hg19_hs37d5_GRCh37", "mm10_GRCm38"]), SHARED_PATH),


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
        gtf=GTFFILE,
    output:
        join("rsemref","{genome}.transcripts.ump"),
    params:
        rname='bl:rsem',
        genome=GENOME,
        prefix=join("rsemref", GENOME),
    threads: int(allocated("threads", "rsem", cluster)),
    container: config['images']['rsem']
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
        "geneinfo.bed",
    params:
        rname='bl:annotate',
        get_gene=join(SCRIPTSDIR, "get_gene_annotate.py"),
        get_isoform=join(SCRIPTSDIR, "get_isoform_annotate.py"),
        make_refFlat=join(SCRIPTSDIR, "make_refFlat.py"),
        make_geneinfo=join(SCRIPTSDIR, "make_geneinfo.py"),
    container: config['images']['build_rnaseq']
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
        gtf=GTFFILE,
    output:
        join("STAR", "2.7.6a", "genes-{readlength}", "SA"),
    params:
        rname='bl:star_rl',
        tmpdir=tmpdir,
    threads: int(allocated("threads", "star_rl", cluster)),
    container: config['images']['arriba']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Create Index for read length
    rl=$(({wildcards.readlength}-1))

    STAR \\
        --runThreadN {threads} \\
        --runMode genomeGenerate \\
        --genomeDir STAR/2.7.6a/genes-{wildcards.readlength} \\
        --genomeFastaFiles {input.fa} \\
        --sjdbGTFfile {input.gtf} \\
        --sjdbOverhang $rl \\
        --outFileNamePrefix STAR/2.7.6a/build_{wildcards.readlength}_ \\
        --outTmpDir ${{tmp}}/tmp_{wildcards.readlength}
    """


if SMALL_GENOME == "True":
    # Build a index that is optimized for 
    # small genomes. For small genomes, the 
    # parameter --genomeSAindexNbases must 
    # to be scaled down, with a typical value 
    # of min(14, log2(GenomeLength)/2 - 1). 
    # For example, for 1 megaBase genome, this 
    # is equal to 9, for 100 kiloBase genome, 
    # this is equal to 7. Using this guidance
    # from the author of STAR, we will dynamically
    # determine what the optimal value should be 
    # based on the provided reference genome size.
    rule star_genome:
        """
        Builds STAR Index to align reads against reference genome without the
        GTF or readlength provided. This index only contain information pertaining
        to the asssembly or reference genome in FASTA format. This indice represents a
        base index from which processing annotations from a GTF and insert junctions
        on the fly. This has the advantage of saving diskspace as an index will not be
        created for a list of predefined readlengths. This rule replaces star_rl above.
        This rule dynamically determine the optimal --genomeSAindexNbases value before 
        running STAR generateGenome based on the size of the provided reference. This is 
        needed for very small reference genomes. 
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
            gtf=GTFFILE,
        output:
            join("STAR", "2.7.6a", "genome", "SA"),
        params:
            rname='bl:star_genome',
            tmpdir=tmpdir,
        threads: int(allocated("threads", "star_genome", cluster)),
        container: config['images']['arriba']
        shell: """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        # Build an index optimized
        # for small reference genomes
        genomeSize=$(grep -v '^>' {input.fa} | awk '{{sum+=length($0)}}; END {{print sum}}')
        # Calculate min(14, log2(GenomeSize)/2 - 1)
        genomeSAindexNbases=$(python -c "import math; print(min(14,int((math.log($genomeSize,2)/2)-1)))")
        echo "Building index with --genomeSAindexNbases $genomeSAindexNbases"
        STAR \\
            --runThreadN {threads} \\
            --runMode genomeGenerate \\
            --genomeSAindexNbases ${{genomeSAindexNbases}} \\
            --genomeDir STAR/2.7.6a/genome \\
            --genomeFastaFiles {input.fa} \\
            --outFileNamePrefix STAR/2.7.6a/build_genome_ \\
            --outTmpDir ${{tmp}}/tmp_genome
        """
else:
    # Build a normal genomic index
    # built with the assumption that
    # the user provided a large genome,
    # i.e. most mammalian genomes like
    # mouse or human
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
            gtf=GTFFILE,
        output:
            join("STAR", "2.7.6a", "genome", "SA")
        params:
            rname='bl:star_genome',
            tmpdir=tmpdir,
        threads: int(allocated("threads", "star_genome", cluster)),
        container: config['images']['arriba']
        shell: """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        STAR \\
            --runThreadN {threads} \\
            --runMode genomeGenerate \\
            --genomeDir STAR/2.7.6a/genome \\
            --genomeFastaFiles {input.fa} \\
            --outFileNamePrefix STAR/2.7.6a/build_genome_ \\
            --outTmpDir ${{tmp}}/tmp_genome
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
        "{genome}.rRNA_interval_list",
    params:
        rname='bl:rRNA_list',
        genome=GENOME,
        create_rRNA=join(SCRIPTSDIR, "create_rRNA_intervals.py"),
    container: config['images']['build_rnaseq']
    shell: """
    python3 {params.create_rRNA} \\
        {input.fa} \\
        {input.gtf} \\
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
        gtf=GTFFILE,
    output:
        "karyoplot_gene_coordinates.txt",
    params:
        rname='bl:karyo_coord',
        get_karyoplot=join(SCRIPTSDIR, "get_karyoplot_gene_coordinates.py"),
    container: config['images']['build_rnaseq']
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
        gtf=GTFFILE,
    output:
        join("karyobeds", "karyobed.bed"),
    params:
        rname='bl:karyo_bed',
        get_karyoplot=join(SCRIPTSDIR, "get_karyoplot_beds.py"),
    container: config['images']['build_rnaseq']
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
        gtf=GTFFILE,
    output:
        "transcripts.protein_coding_only.bed12",
    params:
        rname='bl:tin_ref',
        gtf2protein=join(SCRIPTSDIR, "gtf2protein_coding_genes.py"),
        gene2transcripts=join(SCRIPTSDIR, "gene2transcripts_add_length.py"),
    container: config['images']['build_rnaseq']
    shell: """
    python {params.gtf2protein} {input.gtf} > protein_coding_genes.lst
    gtfToGenePred -ignoreGroupsWithoutExons -genePredExt {input.gtf} genes.gtf.genePred
    genePredToBed genes.gtf.genePred genes.gtf.genePred.bed

    awk -F '\\t' -v OFS='\\t' '{{print $12,$1}}' genes.gtf.genePred \\
        | sort -k1,1n > gene2transcripts

    while read gene; do
        grep "${{gene}}" gene2transcripts;
    done < protein_coding_genes.lst > gene2transcripts.protein_coding_only

    python {params.gene2transcripts} \\
        gene2transcripts.protein_coding_only \\
        genes.gtf.genePred.bed \\
        > gene2transcripts.protein_coding_only.with_len

    sort -k1,1 -k3,3nr gene2transcripts.protein_coding_only.with_len | \\
        awk -F '\\t' '{{if (!seen[$1]) {{seen[$1]++; print $2}}}}' > protein_coding_only.txt

    while read transcript; do
        grep -m1 "${{transcript}}" genes.gtf.genePred.bed;
    done < <(awk -F '.' '{{print $1}}' protein_coding_only.txt) > transcripts.protein_coding_only.bed12

    rm -f "protein_coding_genes.lst" "genes.gtf.genePred" "genes.gtf.genePred.bed" \\
        "gene2transcripts" "gene2transcripts.protein_coding_only" "protein_coding_only.txt" \\
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
        gtf=GTFFILE,
    output:
        "qualimap_info.txt",
    params:
        rname='bl:qualimapinfo',
        generate_qualimap=join(SCRIPTSDIR, "generate_qualimap_ref.py"),
    container: config['images']['build_rnaseq']
    shell: """
    python3 {params.generate_qualimap} \\
        -g {input.gtf} \\
        -f {input.fa} \\
        -o {output} \\
        --ignore-strange-chrom 2> qualimap_error.log
    """


rule fqscreen_db1:
    """
    Downloads fastq screen bowtie2 databases from the OpenOmics public resource bundle.
    Currently, there are fastq screen indices for the following organisms: hg19, mm9, 
    bateria, fungi, virus, univec vector sequences, and rRNA.

    @Input:
        Genomic FASTA file
        Annotation file in GTF format
    @Output:
        Downloaded bowtie2 indices for hg19, mm9, virus, univec, and rRNA
    """
    input:
        fa=REFFA,
        gtf=GTFFILE,
    output:
        fwd=expand(
            join(SHARED_PATH, "fastq_screen_db", "{{genome}}", "{{ref}}.{ext}.bt2"),
            ext=["1", "2", "3", "4"]
        ),
        rev=expand(
            join(SHARED_PATH, "fastq_screen_db", "{{genome}}", "{{ref}}.rev.{ext}.bt2"),
            ext=["1", "2"]
        ),
    params:
        rname='bl:fqscreen_db1',
        tarname=join("fastq_screen_db.{genome}.tar"),
        tarfile=join(SHARED_PATH, "fastq_screen_db.{genome}.tar"),
        outdir=SHARED_PATH,
    container: config['images']['build_rnaseq']
    shell: """
    wget https://hpc.nih.gov/~OpenOmics/common/{params.tarname} -O {params.tarfile}
    tar vxf {params.tarfile} -C {params.outdir} && rm {params.tarfile}
    """


rule fqscreen_db2:
    """
    Downloads fastq screen bowtie2 databases from the OpenOmics public resource bundle.
    Currently, there are fastq screen indices for the following organisms: hg19, mm9, 
    bateria, fungi, virus, univec vector sequences, and rRNA.

    @Input:
        Genomic FASTA file
        Annotation file in GTF format
    @Output:
        Downloaded bowtie2 indices for bateria, fungi
    """
    input:
        fa=REFFA,
        gtf=GTFFILE,
    output:
        fwd=expand(
            join(SHARED_PATH, "fastq_screen_db", "{{genome}}", "{{ref}}.{ext}.bt2l"),
            ext=["1", "2", "3", "4"]
        ),
        rev=expand(
            join(SHARED_PATH, "fastq_screen_db", "{{genome}}", "{{ref}}.rev.{ext}.bt2l"),
            ext=["1", "2"]
        ),
    params:
        rname='bl:fqscreen_db2',
        tarname=join("fastq_screen_db.{genome}.tar"),
        tarfile=join(SHARED_PATH, "fastq_screen_db.{genome}.tar"),
        outdir=SHARED_PATH,
    container: config['images']['build_rnaseq']
    shell: """
    wget https://hpc.nih.gov/~OpenOmics/common/{params.tarname} -O {params.tarfile}
    tar vxf {params.tarfile} -C {params.outdir} && rm {params.tarfile}
    """


rule fqscreen_conf:
    """
    Downloads fastq screen confs from the OpenOmics public resource bundle.
    @Input:
        Genomic FASTA file
        Annotation file in GTF format
    @Output:
        Downloaded bowtie2 indices for bateria, fungi
    """
    input:
        fa=REFFA,
        gtf=GTFFILE,
    output:
        conf1=join(SHARED_PATH, "fastq_screen_db", "fastq_screen_p1.conf"),
        conf2=join(SHARED_PATH, "fastq_screen_db", "fastq_screen_p2.conf"),
    params:
        rname='bl:fqscreen_conf',
        conf1="fastq_screen_p1.conf",
        conf2="fastq_screen_p2.conf",
        new=join(SHARED_PATH).rstrip('/'),
        outdir=join(SHARED_PATH, "fastq_screen_db")
    container: config['images']['build_rnaseq']
    shell: """
    mkdir -p "{params.outdir}"
    wget https://hpc.nih.gov/~OpenOmics/common/{params.conf1} -O {output.conf1}
    wget https://hpc.nih.gov/~OpenOmics/common/{params.conf2} -O {output.conf2}
    sed -i 's@/data/OpenOmics/references/common@{params.new}@g' {output.conf1}
    sed -i 's@/data/OpenOmics/references/common@{params.new}@g' {output.conf2}
    """


rule kraken2_db:
    """
    Downloads kraken2 database from the OpenOmics public resource bundle.
    @Input:
        Genomic FASTA file
        Annotation file in GTF format
    @Output:
        Downloaded kraken2 db for contamination screening
    """
    input:
        fa=REFFA,
        gtf=GTFFILE,
    output:
        expand(join(SHARED_PATH, "k2_pluspf_20241228", "{ref}.k2d"), ref=["hash", "opts", "taxo"]),
    params:
        rname='bl:kraken2_db',
        outfh=join(SHARED_PATH, "k2_pluspf_20241228.tar.gz"),
        tarfile="k2_pluspf_20241228.tar.gz",
        outdir=SHARED_PATH,
    container: config['images']['build_rnaseq']
    shell: """
    wget https://hpc.nih.gov/~OpenOmics/common/{params.tarfile} -O {params.outfh}
    tar zvxf {params.outfh} -C {params.outdir} && rm {params.outfh}
    """


rule arriba_references:
    """
    Downloads arriba references for hg19, hg38, and mm10 from the OpenOmics public resource bundle.
    These references include a blacklist, cytobands, and protein domains for each reference genome.
    @Input:
        Genomic FASTA file
        Annotation file in GTF format
    @Output:
        Downloaded arriba reference files for gene-fusion calling
    """
    input:
        fa=REFFA,
        gtf=GTFFILE,
    output:
        expand(
            join(SHARED_PATH, "arriba", "blacklist_{ref}_v2.0.0.tsv.gz"), 
            ref=["hg38_GRCh38", "hg19_hs37d5_GRCh37", "mm10_GRCm38"]
        ),
    params:
        rname='bl:arriba_refs',
        outfh=join(SHARED_PATH, "arriba_references.tar.gz"),
        tarfile="arriba_references.tar.gz",
        outdir=SHARED_PATH,
    container: config['images']['build_rnaseq']
    shell: """
    wget https://hpc.nih.gov/~OpenOmics/common/{params.tarfile} -O {params.outfh}
    tar zvxf {params.outfh} -C {params.outdir} && rm {params.outfh}
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
        gtf=GTFFILE,
    output:
        json="{genome}.json",
    params:
        rname='bl:jsonmaker',
        workdir=OUTDIR,
        genome=GENOME,
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
            "/data/OpenOmics/references/common/arriba/blacklist_hg19_hs37d5_GRCh37_v2.0.0.tsv.gz"
            refdict["references"]["rnaseq"]["FUSIONCYTOBAND"] = \
            "/data/OpenOmics/references/common/arriba/cytobands_hg19_hs37d5_GRCh37_v2.0.0.tsv"
            refdict["references"]["rnaseq"]["FUSIONPROTDOMAIN"] = \
            "/data/OpenOmics/references/common/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.0.0.gff3"
        elif 'hg38' in params.genome.lower() or \
        'hs38d' in params.genome.lower() or \
        'grch38' in params.genome.lower():
            refdict["references"]["rnaseq"]["FUSIONBLACKLIST"] = \
            "/data/OpenOmics/references/common/arriba/blacklist_hg38_GRCh38_v2.0.0.tsv.gz"
            refdict["references"]["rnaseq"]["FUSIONCYTOBAND"] = \
            "/data/OpenOmics/references/common/arriba/cytobands_hg38_GRCh38_v2.0.0.tsv"
            refdict["references"]["rnaseq"]["FUSIONPROTDOMAIN"] = \
            "/data/OpenOmics/references/common/arriba/protein_domains_hg38_GRCh38_v2.0.0.gff3"
        elif 'mm10' in params.genome.lower() or \
        'grcm38' in params.genome.lower():
            refdict["references"]["rnaseq"]["FUSIONBLACKLIST"] = \
            "/data/OpenOmics/references/common/arriba/blacklist_mm10_GRCm38_v2.0.0.tsv.gz"
            refdict["references"]["rnaseq"]["FUSIONCYTOBAND"] = \
            "/data/OpenOmics/references/common/arriba/cytobands_mm10_GRCm38_v2.0.0.tsv"
            refdict["references"]["rnaseq"]["FUSIONPROTDOMAIN"] = \
            "/data/OpenOmics/references/common/arriba/protein_domains_mm10_GRCm38_v2.0.0.gff3"

        with open(output.json, 'w') as fp:
            json.dump(refdict, fp, indent=4)

