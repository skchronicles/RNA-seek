from snakemake.utils import validate
from scripts.common import abstract_location

# This container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


# Rules common to RNA-seq pipeline, irrespective if the data is single-end or paired-end
rule picard:
    input:
        file1=join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam"),
    output:
        outstar2=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
        outstar2b=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bai"),
        outstar3=join(workpath,log_dir,"{name}.star.duplic")
    params:
        rname='pl:picard',
        sampleName="{name}",
    threads: 6
    envmodules: config['bin'][pfamily]['tool_versions']['PICARDVER'],
    container: "docker://nciccbr/ccbr_picard:v0.0.1"
    shell: """
    java -Xmx110g  -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar AddOrReplaceReadGroups \
        I={input.file1} O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam \
        TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample;
    java -Xmx110g -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates \
        I=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam \
        O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam \
        TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam {output.outstar2};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bai {output.outstar2b};
    sed -i 's/MarkDuplicates/picard.sam.MarkDuplicates/g' {output.outstar3};
    """


rule preseq:
    input:
        bam = join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        ccurve = join(workpath,preseq_dir,"{name}.ccurve"),
    params:
        rname = "pl:preseq",
    envmodules: config['bin'][pfamily]['tool_versions']['PRESEQVER'],
    container: "docker://nciccbr/ccbr_preseq:v0.0.1"
    shell:"""
    preseq c_curve -B -o {output.ccurve} {input.bam}
    """


rule stats:
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
    envmodules:
        config['bin'][pfamily]['tool_versions']['PICARDVER'],
        config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        config['bin'][pfamily]['tool_versions']['PYTHONVER']
    container: "docker://nciccbr/ccbr_rstat:v0.0.1"
    shell: """
    java -Xmx110g -jar ${{PICARDJARPATH}}/picard.jar CollectRnaSeqMetrics \
        REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} \
        RIBOSOMAL_INTERVALS={params.rrnalist} \
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
        TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
    sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
    samtools flagstat {input.file1} > {output.outstar2};
    python3 {params.statscript} {input.file1} >> {output.outstar2}
    """


rule rsem_merge:
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
    container: "docker://nciccbr/ccbr_python:v0.0.1"
    shell: """
    python {params.pythonscript} {params.annotate} {params.inputdir} {params.inputdir}
    sed 's/\\t/|/1' {output.gene_counts_matrix} | \
        sed '1 s/^gene_id|GeneName/symbol/' > {output.reformatted}
    """


rule rseqc:
    input:
        file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        out1=join(workpath,rseqc_dir,"{name}.strand.info"),
        out2=join(workpath,rseqc_dir,"{name}.Rdist.info")
    params:
        bedref=config['references'][pfamily]['BEDREF'],
        rname="pl:rseqc",
    envmodules: config['bin'][pfamily]['tool_versions']['RSEQCVER'],
    container: "docker://nciccbr/ccbr_rseqc_3.0.0:v032219"
    shell: """
    infer_experiment.py -r {params.bedref} -i {input.file1} -s 1000000 > {output.out1}
    read_distribution.py -i {input.file1} -r {params.bedref} > {output.out2}
    """


rule rnaseq_multiqc:
    input:
        expand(join(workpath,rseqc_dir,"{name}.Rdist.info"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,log_dir,"{name}.flagstat.concord.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.star.duplic"),name=samples),
        expand(join(workpath,preseq_dir,"{name}.ccurve"),name=samples),
        expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"),name=samples),
        expand(join(workpath,rseqc_dir,"{name}.Rdist.info"),name=samples),
        qcconfig=abstract_location(config['bin'][pfamily]['CONFMULTIQC']),
    output:
        join(workpath,"Reports","multiqc_report.html"),
        join(workpath,"Reports", "multiqc_matrix.tsv"),
    params:
        rname="pl:multiqc",
        workdir=join(workpath),
        outdir=join(workpath,"Reports"),
        logfiles=join(workpath,"Reports","multiqc_data","*.txt"),
        pyparser=join("workflow", "scripts", "pyparser.py"),
    threads: 2
    envmodules: config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
    container: "docker://nciccbr/ccbr_multiqc_1.9:v0.0.1"
    shell: """
    multiqc -f -c {input.qcconfig} --interactive --outdir {params.outdir} {params.workdir}
    python3 {params.pyparser} {params.logfiles} {params.outdir}
    """
