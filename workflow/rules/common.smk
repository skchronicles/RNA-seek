from snakemake.utils import validate

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

# Helper functions
def check_existence(filename):
    """Checks if file exists on filesystem
    :param filename <str>: Name of file to check
    """
    if not os.path.exists(filename):
        sys.exit("File: {} does not exists!".format(filename))


def check_readaccess(filename):
    """Checks permissions to see if user can read a file
    :param filename <str>: Name of file to check
    """
    check_existence(filename)
    if not os.access(filename,os.R_OK):
        sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))


def check_writeaccess(filename):
    """Checks permissions to see if user can write to a file
    :param filename <str>: Name of file to check
    """
    check_existence(filename)
    if not os.access(filename,os.W_OK):
        sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))


# Rules common to RNA-seq pipeline, irrespective if the data is single-end or paired-end
rule samplecondition:
    input:
        files=expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"), name=samples)
    output:
        out1=join(workpath,star_dir,"sampletable.txt")
    params:
        rname='pl:samplecondition',
        pathprefix=join(workpath,star_dir),
        groups=config['project']['groups']['rgroups'],
        labels=config['project']['groups']['rlabels'],
        allsamples=config['project']['groups']['rsamps'],
        gtffile=config['references'][pfamily]['GTFFILE']
    run:
        with open(output.out1, "w") as out:
            out.write("sampleName\tfileName\tcondition\tlabel\n")
            i=0
            for f in input.files:
                out.write("{}\t".format(params.allsamples[i]))
                out.write("{}/{}.star.count.txt\t".format(params.pathprefix, params.allsamples[i]))
                out.write("{}\t".format(params.groups[i]))
                out.write("{}\n".format(params.labels[i]))
                i=i+1
            out.close()


rule get_strandness:
    input:
        groupsfile=join(workpath,"groups.tab"),
        files=expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
    output:
        outfile=join(workpath,log_dir,"strandness.txt")
    params:
        rname='pl:get_strandness',
        outdir=join(workpath,log_dir),
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        pythonscript=join("workflow", "scripts", "get_strandness.py")
    run:
        import os
        os.chdir(params.outdir)
        check_readaccess(input.groupsfile)
        os.system("module load "+params.pythonver+";python "+params.pythonscript+" "+input.groupsfile+" > "+output.outfile)
        strandfile=open(output.outfile,'r')
        strandness=strandfile.readline().strip()
        strandfile.close()
        A=open(join(workpath,"run.json"),'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['STRANDED']=strandness
        with open(join(workpath,'run.json'),'w') as F:
            json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()


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
    java -Xmx110g  -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar AddOrReplaceReadGroups I={input.file1} O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample;
    java -Xmx110g -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates I=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam {output.outstar2};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bai {output.outstar2b};
    sed -i 's/MarkDuplicates/picard.sam.MarkDuplicates/g' {output.outstar3};
    """


rule preseq:
    params:
        rname = "pl:preseq",
        preseqver=config['bin'][pfamily]['tool_versions']['PRESEQVER'],
    input:
        bam = join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        ccurve = join(workpath,preseq_dir,"{name}.ccurve"),
    shell:"""
    module load {params.preseqver};
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
        picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
        samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        refflat=config['references'][pfamily]['REFFLAT'],
        rrnalist=config['references'][pfamily]['RRNALIST'],
        picardstrand=config['bin'][pfamily]['PICARDSTRAND'],
        statscript=join("workflow", "scripts", "bam_count_concord_stats.py")
    shell: """
    module load R/3.5;
    module load {params.picardver};
    java -Xmx110g -jar $PICARDJARPATH/picard.jar CollectRnaSeqMetrics REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist}  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
    sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
    module load {params.samtoolsver};
    samtools flagstat {input.file1} > {output.outstar2};
    module load python/3.5;
    python {params.statscript} {input.file1} >> {output.outstar2}
    """


rule rnaseq_multiqc:
    input:
        expand(join(workpath,rseqc_dir,"{name}.Rdist.info"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,log_dir,"{name}.flagstat.concord.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.star.duplic"),name=samples),
        expand(join(workpath,preseq_dir,"{name}.ccurve"),name=samples),
        expand(join(workpath,degall_dir,"PcaReport_{dtype}.html"),dtype=dtypes),
    output:
        join(workpath,"Reports","multiqc_report.html")
    params:
        rname="pl:multiqc",
        logsdir=join(workpath,log_dir),
        outdir=join(workpath,"Reports"),
        multiqcver=config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
        qcconfig=config['bin'][pfamily]['CONFMULTIQC']
    threads: 1
    shell: """
    module load {params.multiqcver}
    cd {params.outdir}
    multiqc -f -c {params.qcconfig} --interactive -e cutadapt -d ../
    cd {workpath}/slurmfiles
    multiqc -f --interactive .
    """


rule rsem_merge:
    input:
        files=expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"), name=samples),
        files2=expand(join(workpath,degall_dir,"{name}.RSEM.isoforms.results"), name=samples),
    output:
        join(workpath,degall_dir,"RSEM.genes.FPKM.all_samples.txt"),
        join(workpath,degall_dir,"RSEM.isoforms.FPKM.all_samples.txt"),
    params:
        rname='pl:rsem_merge',
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join("workflow", "scripts", "merge_rsem_results.py"),
    shell: """
    module load {params.pythonver}
    python {params.pythonscript} {params.annotate} {degall_dir} {degall_dir}
    """


rule rsemcounts:
    input:
        files=expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"), name=samples),
        sampletable=join(workpath,star_dir,"sampletable.txt")
    output:
        join(workpath,degall_dir,"RawCountFile_RSEM_genes_filtered.txt"),
    params:
        rname='pl:rsemcounts',
        outdir=join(workpath,degall_dir),
        annotate=config['references'][pfamily]['ANNOTATE'],
        rver=config['bin'][pfamily]['tool_versions']['RVER'],
        rscript=join("workflow", "scripts", "rsemcounts.R")
    shell: """
    cd {params.outdir}
    module load {params.rver}
    Rscript {params.rscript} '{params.outdir}' '{input.files}' '{params.annotate}' '{input.sampletable}'
    """


rule qualicounts:
    input:
        countsmatrix=join(workpath,degall_dir,"RawCountFile_RSEM_genes_filtered.txt"),
        groupsfile=join(workpath,"groups.tab"),
    output:
        outcounts=join(workpath,"QualiMap","RawCountFile_RSEM_genes_filtered_qualimap.txt"),
        globalreport=join(workpath,"QualiMap","GlobalReport.html"),
    params:
        rname='pl:qualicounts',
        sampletable=join(workpath,"QualiMap","qualimap_sample_table.txt"),
        outdir=join(workpath,"QualiMap"),
        info=config['references'][pfamily]['QUALIMAP_INFO'],
    shell: """
    module load qualimap/2.2.1
    # Remove gene symbols from count matrix
    sed 's/|[a-zA-Z0-9]\+//g' {input.countsmatrix} | tail -n +2 > {output.outcounts}
    sed '/^$/d' {input.groupsfile} | awk -v OFS='\\t' '{{print $1, $2,"{output.outcounts}", NR+1}}' > {params.sampletable}
    qualimap counts -d {params.sampletable} -i {params.info} -outdir {params.outdir}
    """


rule rseqc:
    input:
        file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        out1=join(workpath,rseqc_dir,"{name}.strand.info"),
        out4=join(workpath,rseqc_dir,"{name}.Rdist.info")
    params:
        bedref=config['references'][pfamily]['BEDREF'],
        rseqcver=config['bin'][pfamily]['tool_versions']['RSEQCVER'],
        rname="pl:rseqc"
    shell: """
    module load {params.rseqcver}
    cd {rseqc_dir}
    infer_experiment.py -r {params.bedref} -i {input.file1} -s 1000000 > {output.out1}
    read_distribution.py -i {input.file1} -r {params.bedref} > {output.out4}
    """


rule pca:
    input:
        file1=join(workpath,star_dir,"sampletable.txt"),
        file2=join(workpath,degall_dir,"RawCountFile_{dtype}_filtered.txt"),
    output:
        outhtml=join(workpath,degall_dir,"PcaReport_{dtype}.html")
    params:
        rname='pl:pca',
        outdir=join(workpath,degall_dir),
        dtype="{dtype}",
        projectId=config['project']['id'],
        projDesc=config['project']['description'].rstrip('\n'),
        rver=config['bin'][pfamily]['tool_versions']['RVER'],
        scripts_dir=join("workflow", "scripts"),
        rscript1=join("workflow", "scripts", "pcacall.R"),
        rscript2=join("workflow", "scripts", "PcaReport.Rmd"),
    shell: """
    cd {params.outdir}

    module load {params.rver}
    Rscript {params.rscript1} '{params.outdir}' '{output.outhtml}' '{input.file1}' '{input.file2}' '{params.projectId}' '{params.projDesc}' '{params.rscript2}'
    """
