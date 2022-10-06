# Single-end snakemake rules imported in the main Snakefile.

# Pre Alignment Rules
rule validator:
    """
    Quality-control step to ensure the input FastQC files are not corrupted or
    incomplete prior to running the entire workflow. Please note this rule will
    only run if the --use-singularity flag is provided to snakemake.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Log file containing any warnings or errors on file
    """
    input:
        R1=join(workpath,"{name}.R1.fastq.gz"),
    output:
        out1=join(workpath,"rawQC","{name}.validated.R1.fastq.log"),
    priority: 2
    params:
        rname='pl:validator',
        outdir=join(workpath,"rawQC"),
    container: config['images']['fastqvalidator']
    shell: """
    mkdir -p {params.outdir}
    fastQValidator --noeof --file {input.R1} > {output.out1}
    """


rule rawfastqc:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        List of Raw FastQ files (gather)
    @Output:
        List of FastQC reports and zip file containing data quality information
    """
    input:
        R1=join(workpath,"{name}.R1.fastq.gz"),
    output:
        join(workpath,"rawQC","{name}.R1_fastqc.zip"),
    priority: 2
    params:
        rname='pl:rawfastqc',
        outdir=join(workpath,"rawQC"),
    threads: 32
    envmodules: config['bin'][pfamily]['tool_versions']['FASTQCVER']
    container: config['images']['fastqc']
    shell: """
    fastqc {input.R1} -t {threads} -o {params.outdir};
    """

if config['options']['small_rna']:
    # Run STAR with ENCODE's recommendations for small RNA sequencing.
    # Set the min read legth to 
    rule trim_se:
        """
        Data-processing step to remove adapter sequences and perform quality trimming
        prior to alignment the reference genome.  Adapters are composed of synthetic
        sequences and should be removed prior to alignment.
        @Input:
            Raw FastQ file (scatter)
        @Output:
            Trimmed FastQ file
        """
        input:
            infq=join(workpath,"{name}.R1.fastq.gz"),
        output:
            outfq=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"))
        params:
            rname='pl:trim_se',
            # Exposed Parameters: modify config/templates/tools.json to change defaults
            fastawithadaptersetd=config['bin'][pfamily]['tool_parameters']['FASTAWITHADAPTERSETD'],
            leadingquality=config['bin'][pfamily]['tool_parameters']['LEADINGQUALITY'],
            trailingquality=config['bin'][pfamily]['tool_parameters']['TRAILINGQUALITY'],
            minlen=config['bin'][pfamily]['tool_parameters']['MINLEN'],
        threads:32
        envmodules: config['bin'][pfamily]['tool_versions']['CUTADAPTVER']
        container: config['images']['cutadapt']
        shell: """
        cutadapt --nextseq-trim=2 --trim-n \
            -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} \
            -m 16 -b file:{params.fastawithadaptersetd} -j {threads} \
            -o {output.outfq} {input.infq}
        """
else:
    # Use default trimming rule for long RNAs
    rule trim_se:
        """
        Data-processing step to remove adapter sequences and perform quality trimming
        prior to alignment the reference genome.  Adapters are composed of synthetic
        sequences and should be removed prior to alignment.
        @Input:
            Raw FastQ file (scatter)
        @Output:
            Trimmed FastQ file
        """
        input:
            infq=join(workpath,"{name}.R1.fastq.gz"),
        output:
            outfq=temp(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"))
        params:
            rname='pl:trim_se',
            # Exposed Parameters: modify config/templates/tools.json to change defaults
            fastawithadaptersetd=config['bin'][pfamily]['tool_parameters']['FASTAWITHADAPTERSETD'],
            leadingquality=config['bin'][pfamily]['tool_parameters']['LEADINGQUALITY'],
            trailingquality=config['bin'][pfamily]['tool_parameters']['TRAILINGQUALITY'],
            minlen=config['bin'][pfamily]['tool_parameters']['MINLEN'],
        threads:32
        envmodules: config['bin'][pfamily]['tool_versions']['CUTADAPTVER']
        container: config['images']['cutadapt']
        shell: """
        cutadapt --nextseq-trim=2 --trim-n \
            -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} \
            -m {params.minlen} -b file:{params.fastawithadaptersetd} -j {threads} \
            -o {output.outfq} {input.infq}
        """


rule fastqc:
    """
    Quality-control step to assess sequencing quality of the raw data after removing
    adapter sequences. This step is run after trim_pe rule. FastQC is run after adapter
    trimming to evalute if the adapter sequences were properly removed.
    @Input:
        List of Trimmed FastQ files (gather)
    @Output:
        List of FastQC reports and zip file containing data quality information
    """
    input:
        R1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
    output:
        join(workpath,"QC","{name}.R1.trim_fastqc.zip"),
    priority: 2
    params:
        rname='pl:fastqc',
        outdir=join(workpath,"QC"),
        getrl=join("workflow", "scripts", "get_read_length.py"),
    threads: 32
    envmodules: config['bin'][pfamily]['tool_versions']['FASTQCVER']
    container: config['images']['fastqc']
    shell: """
    fastqc {input} -t {threads} -o {params.outdir};
    """

rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
    output:
        out1=join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),
        out2=join(workpath,"FQscreen","{name}.R1.trim_screen.png"),
        out3=join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),
        out4=join(workpath,"FQscreen2","{name}.R1.trim_screen.png")
    params:
        rname='pl:fqscreen',
        outdir = join(workpath,"FQscreen"),
        outdir2 = join(workpath,"FQscreen2"),
        # Exposed Parameters: modify resources/fastq_screen{_2}.conf to change defaults
        # locations to bowtie2 indices
        fastq_screen_config=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG'],
        fastq_screen_config2=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG2'],
    threads: 24
    envmodules:
        config['bin'][pfamily]['tool_versions']['FASTQSCREENVER'],
        config['bin'][pfamily]['tool_versions']['PERLVER'],
        config['bin'][pfamily]['tool_versions']['BOWTIE2VER'],
    container: config['images']['fastq_screen']
    shell: """
    fastq_screen --conf {params.fastq_screen_config} --outdir {params.outdir} \
        --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1}

    fastq_screen --conf {params.fastq_screen_config2} --outdir {params.outdir2} \
        --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1}
    """


rule kraken_se:
    """
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interative krona report
    """
    input:
        fq=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
    output:
        krakenout = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.out.txt"),
        krakentaxa = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.taxa.txt"),
        kronahtml = join(workpath,kraken_dir,"{name}.trim.kraken_bacteria.krona.html"),
    params:
        rname='pl:kraken',
        outdir=join(workpath,kraken_dir),
        bacdb=config['bin'][pfamily]['tool_parameters']['KRAKENBACDB'],
    threads: 24
    envmodules:
        config['bin'][pfamily]['tool_versions']['KRAKENVER'],
        config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
    container: config['images']['kraken']
    shell: """
    # Copy kraken2 db to /lscratch or /dev/shm (RAM-disk) to reduce filesytem strain
    cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
    kdb_base=$(basename {params.bacdb})
    kraken2 --db /lscratch/$SLURM_JOBID/${{kdb_base}} \
        --threads {threads} --report {output.krakentaxa} \
        --output {output.krakenout} \
        --gzip-compressed \
        {input.fq}
    # Generate Krona Report
    cut -f2,3 {output.krakenout} | \
        ktImportTaxonomy - -o {output.kronahtml}
    """

if config['options']['star_2_pass_basic']:
    # Run STAR with per-sample 2-pass mapping using '--twopassMode Basic' option
    # STAR will perform the 1st pass mapping, then it will automatically extract
    # splice junctions, insert them into the genome index, and, finally, re-map
    # all reads in the 2nd mapping pass.
    rule star_basic:
        """
        Data processing step to align reads against reference genome using STAR in
        per sample two-pass basic mode. STAR will perform the 1st pass mapping, then
        it will automatically extract splice junctions, insert them into the genome
        index, and, finally, re-map all reads in the 2nd mapping pass. Agian, Splice
        junctions are detected at a per sample level.
        @Input:
            Trimmed FastQ files (scatter)
        @Output:
            Genomic and transcriptomic BAM file
        """
        input:
            file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        output:
            out1=temp(join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam")),
            out2=join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"),
            out3=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
            out4=join(workpath,star_dir,"{name}.p2.SJ.out.tab"),
            out5=join(workpath,log_dir,"{name}.p2.Log.final.out"),
        params:
            rname='pl:star_basic',
            prefix=join(workpath, star_dir, "{name}.p2"),
            best_rl_script=join("workflow", "scripts", "optimal_read_length.py"),
            # Exposed Parameters: modify config/genomes/{genome}.json to change default
            stardir=config['references'][pfamily]['GENOME_STARDIR'],
            gtffile=config['references'][pfamily]['GTFFILE'],
            # Exposed STAR Parameters: modify config/templates/tools.json to change defaults
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
            nbjuncs=config['bin'][pfamily]['NBJUNCS'],
        threads:32
        envmodules: config['bin'][pfamily]['tool_versions']['STARVER']
        container: config['images']['arriba']
        shell: """
        # Optimal readlength for sjdbOverhang = max(ReadLength) - 1 [Default: 100]
        readlength=$(zcat {input.file1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
        echo "sjdbOverhang for STAR: ${{readlength}}"

        STAR --genomeDir {params.stardir} \
            --outFilterIntronMotifs {params.filterintronmotifs} \
            --outSAMstrandField {params.samstrandfield} \
            --outFilterType {params.filtertype} \
            --outFilterMultimapNmax {params.filtermultimapnmax} \
            --alignSJoverhangMin {params.alignsjoverhangmin} \
            --alignSJDBoverhangMin {params.alignsjdboverhangmin} \
            --outFilterMismatchNmax {params.filtermismatchnmax} \
            --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax} \
            --alignIntronMin {params.alignintronmin} \
            --alignIntronMax {params.alignintronmax} \
            --alignMatesGapMax {params.alignmatesgapmax} \
            --clip3pAdapterSeq {params.adapter1} {params.adapter2} \
            --readFilesIn {input.file1} --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --outSAMunmapped {params.outsamunmapped} \
            --outWigType {params.wigtype} \
            --outWigStrand {params.wigstrand} \
            --twopassMode Basic \
            --sjdbGTFfile {params.gtffile} \
            --limitSjdbInsertNsj {params.nbjuncs} \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --alignEndsProtrude 10 ConcordantPair \
            --peOverlapNbasesMin 10 \
            --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp_{wildcards.name} \
            --sjdbOverhang ${{readlength}}

        mv {params.prefix}.Aligned.toTranscriptome.out.bam {workpath}/{bams_dir};
        mv {params.prefix}.Log.final.out {workpath}/{log_dir}
        """
elif config['options']['small_rna']:
    # Run STAR with ENCODE's recommendations for small RNA sequencing.
    rule star_small:
        """
        Data processing step to align reads against reference genome using STAR using
        ENCODE's recommendations for small RNA. 
        Please see this links for more information:
        https://www.encodeproject.org/pipelines/ENCPL337CSA/
        https://github.com/ENCODE-DCC/long-rna-seq-pipeline/tree/master/dnanexus/small-rna
        @Input:
            Trimmed FastQ files (scatter)
        @Output:
            Genomic and transcriptomic BAM file
        """
        input:
            file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        output:
            out1=temp(join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam")),
            out2=join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"),
            out3=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
            out5=join(workpath,log_dir,"{name}.p2.Log.final.out"),
        params:
            rname='pl:star_small',
            prefix=join(workpath, star_dir, "{name}.p2"),
            best_rl_script=join("workflow", "scripts", "optimal_read_length.py"),
            # Exposed Parameters: modify config/genomes/{genome}.json to change default
            stardir=config['references'][pfamily]['GENOME_STARDIR'],
            gtffile=config['references'][pfamily]['GTFFILE'],
            # Exposed STAR Parameters: modify config/templates/tools.json to change defaults
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
            nbjuncs=config['bin'][pfamily]['NBJUNCS'],
        threads:32
        envmodules: config['bin'][pfamily]['tool_versions']['STARVER']
        container: config['images']['arriba']
        shell: """
        # Optimal readlength for sjdbOverhang = max(ReadLength) - 1 [Default: 100]
        readlength=$(zcat {input.file1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
        echo "sjdbOverhang for STAR: ${{readlength}}"

        STAR --genomeDir {params.stardir} \
            --outFilterMultimapNmax {params.filtermultimapnmax} \
            --alignSJDBoverhangMin 1000 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMatchNmin 16 \
            --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax} \
            --alignIntronMax 1 \
            --clip3pAdapterSeq {params.adapter1} {params.adapter2} \
            --readFilesIn {input.file1} --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --outSAMunmapped Within \
            --sjdbGTFfile {params.gtffile} \
            --limitSjdbInsertNsj {params.nbjuncs} \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM Unsorted \
            --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp_{wildcards.name} \
            --sjdbOverhang ${{readlength}}

        # SAMtools sort (uses less memory than STAR SortedByCoordinate)
        samtools sort -@ {threads} \
            -m 2G -T /lscratch/${{SLURM_JOB_ID}}/SORTtmp_{wildcards.name} \
            -O bam {params.prefix}.Aligned.out.bam \
            > {output.out1}

        rm {params.prefix}.Aligned.out.bam
        mv {params.prefix}.Aligned.toTranscriptome.out.bam {workpath}/{bams_dir};
        mv {params.prefix}.Log.final.out {workpath}/{log_dir}
        """
else:
    # Run STAR with multi-sample 2-pass mapping
    # For a study with multiple samples, it is recommended to collect 1st pass
    # splice junctions from all samples and provide them to the second pass of STAR
    rule star1p:
        """
        Data processing step to align reads against reference genome using STAR in
        two-pass mode. STAR is run in a two-pass mode for enhanced detection of reads
        mapping to novel splice junctions. This rule represents the first pass of STAR.
        @Input:
            Trimmed FastQ files (scatter)
        @Output:
            Logfile containing splice-junctions detected by STAR
        """
        input:
            file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        output:
            out1=join(workpath,star_dir,"{name}.SJ.out.tab"),
            out3=temp(join(workpath,star_dir,"{name}.Aligned.out.bam")),
        params:
            rname='pl:star1p',
            prefix=join(workpath, star_dir, "{name}"),
            best_rl_script=join("workflow", "scripts", "optimal_read_length.py"),
            # Exposed Parameters: modify config/genomes/{genome}.json to change default
            stardir=config['references'][pfamily]['GENOME_STARDIR'],
            gtffile=config['references'][pfamily]['GTFFILE'],
            # Exposed STAR Parameters: modify config/templates/tools.json to change defaults
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
        envmodules: config['bin'][pfamily]['tool_versions']['STARVER']
        container: config['images']['arriba']
        shell: """
        # Optimal readlength for sjdbOverhang = max(ReadLength) - 1 [Default: 100]
        readlength=$(zcat {input.file1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
        echo "sjdbOverhang for STAR: ${{readlength}}"

        STAR --genomeDir {params.stardir} \
            --outFilterIntronMotifs {params.filterintronmotifs} \
            --outSAMstrandField {params.samstrandfield}  \
            --outFilterType {params.filtertype} \
            --outFilterMultimapNmax {params.filtermultimapnmax} \
            --alignSJoverhangMin {params.alignsjoverhangmin} \
            --alignSJDBoverhangMin {params.alignsjdboverhangmin} \
            --outFilterMismatchNmax {params.filtermismatchnmax} \
            --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax} \
            --alignIntronMin {params.alignintronmin} \
            --alignIntronMax {params.alignintronmax} \
            --alignMatesGapMax {params.alignmatesgapmax} \
            --clip3pAdapterSeq {params.adapter1} {params.adapter2} \
            --readFilesIn {input.file1} --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --outSAMtype BAM Unsorted \
            --alignEndsProtrude 10 ConcordantPair \
            --peOverlapNbasesMin 10 \
            --sjdbGTFfile {params.gtffile} \
            --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp_{wildcards.name} \
            --sjdbOverhang ${{readlength}}
        """


    rule sjdb:
        """
        Aggregation step to collect the set of all novel junctions that were detected
        in the first-pass of STAR. These splice junctions will be used to re-build the
        genomic indices.
        @Input:
            Logfiles containing splice-junctions (gather)
        @Output:
            Logfile containing the set of all splice junctions across all samples
        """
        input:
            files=expand(join(workpath,star_dir,"{name}.SJ.out.tab"), name=samples),
        output:
            out1=join(workpath,star_dir,"uniq.filtered.SJ.out.tab"),
        params:
            rname='pl:sjdb'
        shell: """
        cat {input.files} | \
            sort | uniq | \
            awk -F '\\t' '{{if ($5>0 && $6==1) {{print}}}}'| \
            cut -f1-4 | sort | uniq | \
        grep \"^chr\"|grep -v \"^chrM\" > {output.out1}
        """


    rule star2p:
        """
        Data processing step to align reads against reference genome using STAR in
        two-pass mode. This step represents the second step of STAR. Here set of splice
        all novel junctions that were detected in the first-pass of STAR are then inserted
        into the genome indices. In this second-pass, all reads are re-mapped using
        the annotated junctions from the GTF file and novel junctions that were
        detected in the first-pass of STAR.
        @Input:
            Trimmed FastQ files (scatter)
        @Output:
            Genomic and transcriptomic BAM file
        """
        input:
            file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
            tab=join(workpath,star_dir,"uniq.filtered.SJ.out.tab"),
        output:
            out1=temp(join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam")),
            out2=join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"),
            out3=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
            out4=join(workpath,star_dir,"{name}.p2.SJ.out.tab"),
            out5=join(workpath,log_dir,"{name}.p2.Log.final.out"),
        params:
            rname='pl:star2p',
            prefix=join(workpath, star_dir, "{name}.p2"),
            best_rl_script=join("workflow", "scripts", "optimal_read_length.py"),
            # Exposed Parameters: modify config/genomes/{genome}.json to change default
            stardir=config['references'][pfamily]['GENOME_STARDIR'],
            gtffile=config['references'][pfamily]['GTFFILE'],
            # Exposed STAR Parameters: modify config/templates/tools.json to change defaults
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
            nbjuncs=config['bin'][pfamily]['NBJUNCS'],
        threads:32
        envmodules: config['bin'][pfamily]['tool_versions']['STARVER']
        container: config['images']['arriba']
        shell: """
        # Optimal readlength for sjdbOverhang = max(ReadLength) - 1 [Default: 100]
        readlength=$(zcat {input.file1} | awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; END {{print maxlen-1}}')
        echo "sjdbOverhang for STAR: ${{readlength}}"

        STAR --genomeDir {params.stardir} \
            --outFilterIntronMotifs {params.filterintronmotifs} \
            --outSAMstrandField {params.samstrandfield} \
            --outFilterType {params.filtertype} \
            --outFilterMultimapNmax {params.filtermultimapnmax} \
            --alignSJoverhangMin {params.alignsjoverhangmin} \
            --alignSJDBoverhangMin {params.alignsjdboverhangmin} \
            --outFilterMismatchNmax {params.filtermismatchnmax} \
            --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax} \
            --alignIntronMin {params.alignintronmin} \
            --alignIntronMax {params.alignintronmax} \
            --alignMatesGapMax {params.alignmatesgapmax} \
            --clip3pAdapterSeq {params.adapter1} {params.adapter2} \
            --readFilesIn {input.file1} --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --outSAMunmapped {params.outsamunmapped} \
            --outWigType {params.wigtype} \
            --outWigStrand {params.wigstrand} \
            --sjdbFileChrStartEnd {input.tab} \
            --sjdbGTFfile {params.gtffile} \
            --limitSjdbInsertNsj {params.nbjuncs} \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --alignEndsProtrude 10 ConcordantPair \
            --peOverlapNbasesMin 10 \
            --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp_{wildcards.name} \
            --sjdbOverhang ${{readlength}}

        mv {params.prefix}.Aligned.toTranscriptome.out.bam {workpath}/{bams_dir};
        mv {params.prefix}.Log.final.out {workpath}/{log_dir}
        """


# Post Alignment Rules
rule rsem:
    """
    Data processing step to quantify gene and isoform counts.
    @Input:
        Transcriptomic BAM file (scatter)
    @Output:
        Gene and isoform counts
    """
    input:
        file1=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        file2=join(workpath,rseqc_dir,"{name}.strand.info")
    output:
        out1=join(workpath,degall_dir,"{name}.RSEM.genes.results"),
        out2=join(workpath,degall_dir,"{name}.RSEM.isoforms.results"),
    params:
        rname='pl:rsem',
        prefix=join(workpath,degall_dir,"{name}.RSEM"),
        rsemref=config['references'][pfamily]['RSEMREF'],
        annotate=config['references'][pfamily]['ANNOTATE'],
    threads: 16
    envmodules:
        config['bin'][pfamily]['tool_versions']['RSEMVER'],
        config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    container: config['images']['rsem']
    shell: """
    # Get strandedness to calculate Forward Probability
    fp=$(tail -n1 {input.file2} | awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}')

    echo "Forward Probability Passed to RSEM: $fp"
    rsem-calculate-expression --no-bam-output --calc-ci --seed 12345 \
        --bam -p {threads}  {input.file1} {params.rsemref} {params.prefix} --time \
        --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
    """


rule bam2bw_rnaseq_se:
    """
    Summarization step to convert a bam file into forward and reverse strand
    bigwig files suitable for a genomic track viewer like IGV.
    @Input:
        Sorted, duplicate marked genomic BAM file (scatter)
    @Output:
        Forward and reverse strand BigWig files
    """
    input:
        bam=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
        strandinfo=join(workpath,rseqc_dir,"{name}.strand.info")
    output:
        fbw=join(workpath,bams_dir,"{name}.fwd.bw"),
        rbw=join(workpath,bams_dir,"{name}.rev.bw")
    params:
        rname='pl:bam2bw',
        outprefix=join(workpath,bams_dir,"{name}"),
        bashscript=join("workflow", "scripts", "bam2strandedbw.se.sh")
    threads: 2
    envmodules:
        config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
        config['bin'][pfamily]['tool_versions']['BEDTOOLSVER'],
        config['bin'][pfamily]['tool_versions']['UCSCVER'],
        config['bin'][pfamily]['tool_versions']['PARALLELVER'],
    container: config['images']['bam2strandedbw']
    shell: """
    bash {params.bashscript} {input.bam} {params.outprefix}

    # reverse files if method is not dUTP/NSR/NNSR ... ie, R1 in the direction of RNA strand.
    strandinfo=`tail -n1 {input.strandinfo} | awk '{{print $NF}}'`
    if [ `echo "$strandinfo < 0.25"|bc` -eq 1 ]; then
        mv {output.fbw} {output.fbw}.tmp
        mv {output.rbw} {output.fbw}
        mv {output.fbw}.tmp {output.rbw}
    fi
    """


rule rnaseq_multiqc:
    """
    Reporting step to aggregate sample statistics and quality-control information
    across all samples. This will be one of the last steps of the pipeline. The inputs
    listed here are to ensure that this step runs last. During runtime, MultiQC will
    recurively crawl through the working directory and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    input:
        expand(join(workpath,rseqc_dir,"{name}.Rdist.info"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,log_dir,"{name}.flagstat.concord.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.star.duplic"),name=samples),
        expand(join(workpath,preseq_dir,"{name}.ccurve"),name=samples),
        expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"),name=samples),
        expand(join(workpath,rseqc_dir,"{name}.Rdist.info"),name=samples),
        expand(join(workpath,rseqc_dir,"{name}.star_rg_added.sorted.dmark.summary.txt"),name=samples),
        expand(join(workpath,"rawQC","{name}.fastq.info.txt"),name=samples),
        fqinfo=expand(join(workpath,"rawQC","{name}.fastq.info.txt"),name=samples),
        tins=expand(join(workpath,rseqc_dir,"{name}.star_rg_added.sorted.dmark.summary.txt"),name=samples)
    output:
        join(workpath,"Reports","multiqc_report.html"),
        join(workpath,"Reports", "multiqc_matrix.tsv"),
        medtins=join(workpath,"Reports","multiqc_data", "rseqc_median_tin.txt"),
        fclanes=join(workpath,"Reports","multiqc_data", "fastq_flowcell_lanes.txt"),
    params:
        rname="pl:multiqc",
        workdir=join(workpath),
        outdir=join(workpath,"Reports"),
        logfiles=join(workpath,"Reports","multiqc_data","*.txt"),
        pyparser=join("workflow", "scripts", "pyparser.py"),
        qcconfig=config['bin'][pfamily]['CONFMULTIQC'],
    threads: 2
    envmodules: config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
    container: config['images']['multiqc']
    shell: """
    multiqc --ignore '*/.singularity/*' -f -c {params.qcconfig} --interactive --outdir {params.outdir} {params.workdir}

    # Parse RSeQC Median TINs
    echo -e "Sample\\tmedian_tin" > {output.medtins}
    find {params.workdir} -name '*.star_rg_added.sorted.dmark.summary.txt' -exec cut -f1,3 {{}} \\; | \
        grep -v '^Bam_file' | \
        awk -F '\\t' '{{printf "%s\\t%.3f\\n", $1,$2}}' >> {output.medtins}

    # Parse Flowcell and Lane information
    echo -e "Sample\\tflowcell_lanes" > {output.fclanes}
    find {params.workdir} -name '*.fastq.info.txt' -exec awk -F '\\t' -v OFS='\\t' 'NR==2 {{print $1,$5}}' {{}} \\; \
        >> {output.fclanes}

    python3 {params.pyparser} {params.logfiles} {params.outdir}
    """
