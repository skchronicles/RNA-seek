from snakemake.utils import validate
from scripts.common import (
    abstract_location,
    allocated
)

# This container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: config['images']['miniconda']


# Rules common to RNA-seq pipeline, irrespective if the data is single-end or paired-end
rule fc_lane:
    """
    Quality-control step to get flowcell and lane information from FastQ file.
    FastQ files generated with older versions of Casava or downloaded from
    SRA have a different format than newer FastQ files generated with the
    current version of Casava. It is worth noting that FastQ files downloaded from SRA
    or FastQ files generated with Casava version < 1.8 do not have Flowcell
    IDs in its sequence indentifer. If a FastQ file does not have Flowcell IDs,
    the Machine or Instrument ID is grabbed instead.
    @Input:
        Raw FastQ R1 file (scatter)
    @Output:
        Text file containing information about the FastQ file
    """
    input:
        R1=join(workpath,"{name}.R1.fastq.gz"),
    output:
        fqinfo=join(workpath,"rawQC","{name}.fastq.info.txt")
    params:
        rname='pl:fc_lane',
        get_flowcell_lanes=join("workflow", "scripts", "get_flowcell_lanes.py"),
    envmodules: config['bin'][pfamily]['tool_versions']['PYTHONVER']
    container: config['images']['python']
    shell: """
    python {params.get_flowcell_lanes} {input.R1} {wildcards.name} > {output.fqinfo}
    """


rule picard:
    """
    Data processing and quality-control step to add read groups and mark duplicate reads.
    @Input:
        Genomic BAM file (scatter)
    @Output:
        Read groups added, duplicate marked genomic BAM file (scatter)
    """
    input:
        file1=join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam"),
    output:
        bam=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
        bai=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam.bai"),
        metrics=join(workpath,log_dir,"{name}.star.duplic")
    params:
        rname='pl:picard',
        sampleName="{name}",
        memory=allocated("mem", "picard", cluster).lower().rstrip('g'),
    threads: int(allocated("threads", "picard", cluster))-1,
    envmodules: config['bin'][pfamily]['tool_versions']['PICARDVER'],
    container: config['images']['picard']
    shell: """
    java -Xmx{params.memory}g  -XX:ParallelGCThreads={threads} -jar ${{PICARDJARPATH}}/picard.jar AddOrReplaceReadGroups \
        I={input.file1} O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam \
        TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample VALIDATION_STRINGENCY=SILENT;
    java -Xmx{params.memory}g -XX:ParallelGCThreads={threads} -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates \
        I=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam \
        O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam \
        TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.metrics};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam {output.bam};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bai {output.bai};
    sed -i 's/MarkDuplicates/picard.sam.MarkDuplicates/g' {output.metrics};
    """


rule preseq:
    """
    Quality step to estimate library complexity. Low library complexity may indicate
    an issue with library preparation or sample storage (FFPE samples) where very
    little input RNA was over-amplified or the sample may be highly degraded.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        Logfile containing library complexity information
    """
    input:
        bam = join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        ccurve = join(workpath,preseq_dir,"{name}.ccurve"),
    params:
        rname = "pl:preseq",
    envmodules: config['bin'][pfamily]['tool_versions']['PRESEQVER'],
    container: config['images']['preseq']
    shell:"""
    preseq c_curve -B -o {output.ccurve} {input.bam}
    """


rule qualibam:
    """
    Quality-control step to assess various post-alignment metrics and a secondary
    method to calculate insert size.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        Report containing post-aligment QC metrics
    """
    input:
        bamfile=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        report=join(workpath,"QualiMap","{name}","qualimapReport.html"),
        results=join(workpath,"QualiMap","{name}","genome_results.txt"),
    params:
        rname='pl:qualibam',
        outdir=join(workpath,"QualiMap","{name}"),
        gtfFile=config['references'][pfamily]['GTFFILE'],
        memory=int(allocated("mem", "qualibam", cluster).lower().rstrip('g'))-1,
    threads: int(allocated("threads", "qualibam", cluster)),
    envmodules: config['bin'][pfamily]['tool_versions']['QUALIMAPVER']
    container: config['images']['qualimap']
    shell: """
    unset DISPLAY;
    qualimap bamqc -bam {input.bamfile} --feature-file {params.gtfFile} \
        -outdir {params.outdir} -nt {threads} --java-mem-size={params.memory}G
    """


rule stats:
    """
    Quality step to run Picard CollectRnaSeqMetrics. This step collects alignment
    metrics which are specific to cDNA fragments originating from RNA. It describes
    the distribution of reads aligning to various sub-features like the percentage
    of reads aligning to UTRs, intronic, coding, and intergenic regions. Other metrics
    include the median coverage (depth), the ratios of 5 prime /3 prime-biases,
    and the numbers of reads with the correct/incorrect strand designation.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        Logfiles containing RNA-specific alignment metrics
    """
    input:
        file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        outstar1=join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),
        outstar2=join(workpath,log_dir,"{name}.flagstat.concord.txt"),
    params:
        rname='pl:stats',
        refflat=config['references'][pfamily]['REFFLAT'],
        rrnalist=config['references'][pfamily]['RRNALIST'],
        picardstrand=config['bin'][pfamily]['PICARDSTRAND'],
        statscript=join("workflow", "scripts", "bam_count_concord_stats.py"),
        memory=allocated("mem", "stats", cluster).lower().rstrip('g'),
    envmodules:
        config['bin'][pfamily]['tool_versions']['PICARDVER'],
        config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        config['bin'][pfamily]['tool_versions']['PYTHONVER']
    container: config['images']['rstat']
    shell: """
    java -Xmx{params.memory}g -jar ${{PICARDJARPATH}}/picard.jar CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} \
        RIBOSOMAL_INTERVALS={params.rrnalist} \
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
        TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
    sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
    samtools flagstat {input.file1} > {output.outstar2};
    python3 {params.statscript} {input.file1} >> {output.outstar2}
    """


rule rsem_merge:
    """Data processing step to merge the gene and isoform counts for each sample
    into count matrices.
    @Input:
        List of per sample gene and isoform counts (gather)
    @Output:
        Expected Count, FPKM, and TPM count matrices for genes and isoforms
    """
    input:
        files=expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"), name=samples),
        files2=expand(join(workpath,degall_dir,"{name}.RSEM.isoforms.results"), name=samples),
    output:
        gene_counts_matrix=join(workpath,degall_dir,"RSEM.genes.expected_count.all_samples.txt"),
        gene_fpkm_matrix=join(workpath,degall_dir,"RSEM.genes.FPKM.all_samples.txt"),
        isoform_fpkm_matrix=join(workpath,degall_dir,"RSEM.isoforms.FPKM.all_samples.txt"),
        reformatted=join(workpath,degall_dir,"RSEM_genes_expected_counts.tsv"),
    params:
        rname='pl:rsem_merge',
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join("workflow", "scripts", "merge_rsem_results.py"),
        inputdir=join(workpath, degall_dir)
    envmodules: config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    container: config['images']['python']
    shell: """
    python {params.pythonscript} {params.annotate} {params.inputdir} {params.inputdir}
    sed 's/\\t/|/1' {output.gene_counts_matrix} | \
        sed '1 s/^gene_id|GeneName/symbol/' > {output.reformatted}
    """


rule rseqc:
    """
    Quality-control step to infer stranded-ness and read distributions over
    specific genomic features.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        RSeQC logfile containing strand information and read distributions
    """
    input:
        file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        out1=join(workpath,rseqc_dir,"{name}.strand.info"),
        out2=join(workpath,rseqc_dir,"{name}.Rdist.info")
    params:
        bedref=config['references'][pfamily]['BEDREF'],
        rname="pl:rseqc",
    envmodules: config['bin'][pfamily]['tool_versions']['RSEQCVER'],
    container: config['images']['rseqc']
    shell: """
    infer_experiment.py -r {params.bedref} -i {input.file1} -s 1000000 > {output.out1}
    read_distribution.py -i {input.file1} -r {params.bedref} > {output.out2}
    """


rule tin:
    """
    Quality-control step to infer RNA integrity at the transcript level.
    TINs (transcript integrity numbers) are calculated for all canoncial
    protein-coding transcripts. TIN is analogous to a computionally derived
    RIN value. From the docs: requires a sort and indexed bam file.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        RSeQC logfiles containing transcript integrity number information
    """
    input:
        bam=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
        bai=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam.bai"),
        control=join(workpath,rseqc_dir,"{name}.Rdist.info")
    output:
        out1=join(workpath,rseqc_dir,"{name}.star_rg_added.sorted.dmark.tin.xls"),
        out2=join(workpath,rseqc_dir,"{name}.star_rg_added.sorted.dmark.summary.txt")
    params:
        bedref=config['references'][pfamily]['TINREF'],
        outdir=join(workpath,rseqc_dir),
        rname="pl:tin",
    envmodules: config['bin'][pfamily]['tool_versions']['RSEQCVER'],
    container: config['images']['rseqc']
    shell: """
    # tin.py writes to current working directory
    cd {params.outdir}
    tin.py -i {input.bam} -r {params.bedref}
    """


rule tin_merge:
    """
    Data processing step to merge the transcript integrity numbers (TINs) for
    each sample into a matrix.
    @Input:
        RSeQC transcript integrity numbers (gather)
    @Output:
        TIN value matrix (combined_TIN.tsv)
    """
    input:
        tins=expand(join(workpath,rseqc_dir,"{name}.star_rg_added.sorted.dmark.tin.xls"), name=samples)
    output:
        matrix=join(workpath,degall_dir,"combined_TIN.tsv")
    params:
        rname="pl:tin_merge",
        create_matrix=join("workflow", "scripts", "create_tin_matrix.py")
    envmodules: config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    container: config['images']['python']
    shell: """
    python {params.create_matrix} {input.tins} > {output.matrix}
    """
