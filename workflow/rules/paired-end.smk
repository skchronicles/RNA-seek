# Paired-end snakemake rules imported in the main Snakefile.

# Pre Alignment Rules
rule rawfastqc:
    input:
        expand(join(workpath,"{name}.R1.fastq.gz"), name=samples),
        expand(join(workpath,"{name}.R2.fastq.gz"), name=samples)
    output:
        join(workpath,"rawQC")
    priority: 2
    params:
        rname='pl:rawfastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
    threads: 32
    shell: """
    mkdir -p {output};
    module load {params.fastqcver};
    fastqc {input} -t {threads} -o {output};
    """


rule trim_pe:
    input:
        file1=join(workpath,"{name}.R1."+config['project']['filetype']),
        file2=join(workpath,"{name}.R2."+config['project']['filetype']),
    output:
        out1=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")),
        out2=temp(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"))
    params:
        rname='pl:trim_pe',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        cutadaptver=config['bin'][pfamily]['tool_versions']['CUTADAPTVER'],
        fastawithadaptersetd=config['bin'][pfamily]['tool_parameters']['FASTAWITHADAPTERSETD'],
        leadingquality=config['bin'][pfamily]['tool_parameters']['LEADINGQUALITY'],
        trailingquality=config['bin'][pfamily]['tool_parameters']['TRAILINGQUALITY'],
        minlen=config['bin'][pfamily]['tool_parameters']['MINLEN'],
    threads:32
    shell: """
    module load {params.cutadaptver};
    sample=`echo {input.file1}|awk -F "/" '{{print $NF}}'|awk -F ".R1.fastq" '{{print $1}}'`
    cutadapt --pair-filter=any --nextseq-trim=2 --trim-n -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} -m {params.minlen}:{params.minlen} -b file:{params.fastawithadaptersetd} -B file:{params.fastawithadaptersetd} -j {threads} -o {output.out1} -p {output.out2} {input.file1} {input.file2}
    """


rule fastqc:
    input:
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"), name=samples),
        expand(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"), name=samples)
    output:
        join(workpath,"QC")
    priority: 2
    params:
        rname='pl:fastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
    threads: 32
    shell: """
    mkdir -p {output};
    module load {params.fastqcver};
    fastqc {input} -t {threads} -o {output};
    module load python/3.5;
    python Scripts/get_read_length.py {output} > {output}/readlength.txt  2> {output}/readlength.err
    """


rule fastq_screen:
    input:
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
    output:
        out1=join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),
        out2=join(workpath,"FQscreen","{name}.R1.trim_screen.png"),
        out3=join(workpath,"FQscreen","{name}.R2.trim_screen.txt"),
        out4=join(workpath,"FQscreen","{name}.R2.trim_screen.png"),
        out5=join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),
        out6=join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),
        out7=join(workpath,"FQscreen2","{name}.R2.trim_screen.txt"),
        out8=join(workpath,"FQscreen2","{name}.R2.trim_screen.png")
    params:
        rname='pl:fqscreen',
        batch='--cpus-per-task=24 --mem=64g --time=10:00:00',
        fastq_screen=config['bin'][pfamily]['tool_versions']['FASTQ_SCREEN'],
        outdir = join(workpath,"FQscreen"),
        outdir2 = join(workpath,"FQscreen2"),
        fastq_screen_config=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG'],
        fastq_screen_config2=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG2'],
        perlver=config['bin'][pfamily]['tool_versions']['PERLVER'],
        bowtie2ver=config['bin'][pfamily]['tool_versions']['BOWTIE2VER'],
    threads: 24
    shell: """
    module load {params.bowtie2ver};
    module load {params.perlver};
    {params.fastq_screen} --conf {params.fastq_screen_config} --outdir {params.outdir} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1} {input.file2}
    {params.fastq_screen} --conf {params.fastq_screen_config2} --outdir {params.outdir2} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1} {input.file2}
    """

rule kraken_pe:
    input:
        fq1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        fq2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
    output:
        krakentaxa = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),
        kronahtml = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),
    params:
        rname='pl:kraken',
        prefix = "{name}",
        outdir=join(workpath,kraken_dir),
        bacdb=config['bin'][pfamily]['tool_parameters']['KRAKENBACDB'],
        krakenver=config['bin'][pfamily]['tool_versions']['KRAKENVER'],
        kronatoolsver=config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
    threads: 24
    shell: """
    module load {params.krakenver};
    module load {params.kronatoolsver};
    if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
    cd /lscratch/$SLURM_JOBID;
    cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
    kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload --paired {input.fq1} {input.fq2}
    kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
    cut -f2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml
    mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
    mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
    """


rule star1p:
    input:
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
        qcdir=join(workpath,"QC"),
    output:
        out1= join(workpath,star_dir,"{name}.SJ.out.tab"),
        out3= temp(join(workpath,star_dir,"{name}.Aligned.out.bam")),
    params:
        rname='pl:star1p',prefix="{name}",
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        starver=config['bin'][pfamily]['tool_versions']['STARVER'],
        stardir=config['references'][pfamily]['STARDIR'],
        filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],
        samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],
        filtertype=config['bin'][pfamily]['FILTERTYPE'],
        filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],
        alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],
        alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],
        filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],
        filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],
        alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],
        alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],
        alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],
        adapter1=config['bin'][pfamily]['ADAPTER1'],
        adapter2=config['bin'][pfamily]['ADAPTER2'],
    threads: 32
    run:
        import glob,json
        rl=int(open(join(input.qcdir,"readlength.txt")).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="cd {star_dir}; module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+" --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMtype BAM Unsorted --alignEndsProtrude 10 ConcordantPair --peOverlapNbasesMin 10"
        shell(cmd)
        A=open(join(workpath,"run.json"),'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['SJDBOVERHANG']=str(bestdbrl)
        config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
        config['project']['READLENGTH']=str(rl+1)
        with open(join(workpath,'run.json'),'w') as F:
          json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()


rule sjdb:
    input:
        files=expand(join(workpath,star_dir,"{name}.SJ.out.tab"), name=samples)
    output:
        out1=join(workpath,star_dir,"uniq.filtered.SJ.out.tab")
    params:
        rname='pl:sjdb'
    shell: """
    cat {input.files} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort|uniq|grep \"^chr\"|grep -v \"^chrM\" > {output.out1}
    """


rule star2p:
    input:
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
        tab=join(workpath,star_dir,"uniq.filtered.SJ.out.tab"),
        qcdir=join(workpath,"QC"),
    output:
        out1=temp(join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam")),
        out2=join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"),
        out3=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        out4=join(workpath,star_dir,"{name}.p2.SJ.out.tab"),
        out5=join(workpath,log_dir,"{name}.p2.Log.final.out"),
    params:
        rname='pl:star2p',
        prefix="{name}.p2",
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        starver=config['bin'][pfamily]['tool_versions']['STARVER'],
        filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],
        samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],
        filtertype=config['bin'][pfamily]['FILTERTYPE'],
        filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],
        alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],
        alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],
        filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],
        filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],
        alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],
        alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],
        alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],
        adapter1=config['bin'][pfamily]['ADAPTER1'],
        adapter2=config['bin'][pfamily]['ADAPTER2'],
        outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],
        wigtype=config['bin'][pfamily]['WIGTYPE'],
        wigstrand=config['bin'][pfamily]['WIGSTRAND'],
        gtffile=config['references'][pfamily]['GTFFILE'],
        nbjuncs=config['bin'][pfamily]['NBJUNCS'],
        stardir=config['references'][pfamily]['STARDIR'],
    threads:32
    run:
        import glob,os
        rl=int(open(join(input.qcdir,"readlength.txt")).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="cd {star_dir}; module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 10 ConcordantPair --peOverlapNbasesMin 10;"
        shell(cmd)
        cmd="sleep 120;cd {workpath};mv {workpath}/{star_dir}/{params.prefix}.Aligned.toTranscriptome.out.bam {workpath}/{bams_dir}; mv {workpath}/{star_dir}/{params.prefix}.Log.final.out {workpath}/{log_dir}"
        shell(cmd)


rule qualibam:
    input:
        bamfile=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    output:
        report=join(workpath,"QualiMap","{name}","qualimapReport.html"),
    params:
        rname='pl:qualibam',
        outdir=join(workpath,"QualiMap","{name}"),
        gtfFile=config['references'][pfamily]['GTFFILE'],
    shell: """
    module load qualimap/2.2.1
    unset DISPLAY;qualimap bamqc -bam {input.bamfile} --feature-file {params.gtfFile} -outdir {params.outdir} -nt $SLURM_CPUS_PER_TASK --java-mem-size=11G
    """

# Post Alignment Rules
rule rsem:
    input:
        file1=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        file2=join(workpath,rseqc_dir,"{name}.strand.info")
    output:
        out1=join(workpath,degall_dir,"{name}.RSEM.genes.results"),
        out2=join(workpath,degall_dir,"{name}.RSEM.isoforms.results"),
    params:
        rname='pl:rsem',
        prefix="{name}.RSEM",
        outdir=join(workpath,degall_dir),
        batch='--cpus-per-task=16 --mem=32g --time=24:00:00',
        rsemref=config['references'][pfamily]['RSEMREF'],
        rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
    threads: 16
    shell: """
    if [ ! -d {params.outdir} ]; then mkdir {params.outdir}; fi
    cd {params.outdir}
    module load {params.rsemver}
    fp=`tail -n1 {input.file2} |awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}'`
    echo $fp
    rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p {threads}  {input.file1} {params.rsemref} {params.prefix} --time --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=$fp --estimate-rspd
    """


rule bam2bw_rnaseq_pe:
    input:
        bam=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
        strandinfo=join(workpath,rseqc_dir,"{name}.strand.info")
    output:
        fbw=join(workpath,bams_dir,"{name}.fwd.bw"),
        rbw=join(workpath,bams_dir,"{name}.rev.bw")
    params:
        rname='pl:bam2bw',
        prefix="{name}",
        bashscript=join(workpath,"Scripts","bam2strandedbw.pe.sh")
    threads: 4
    shell: """
    sh {params.bashscript} {input.bam}

    # reverse files if method is not dUTP/NSR/NNSR ... ie, R1 in the direction of RNA strand.
    strandinfo=`tail -n1 {input.strandinfo}|awk '{{print $NF}}'`
    if [ `echo "$strandinfo < 0.25"|bc` -eq 1 ];then
    mv {output.fbw} {output.fbw}.tmp
    mv {output.rbw} {output.fbw}
    mv {output.fbw}.tmp {output.rbw}
    fi
    """
