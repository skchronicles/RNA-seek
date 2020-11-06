# Group-info snakemake rules imported in the main Snakefile.

# Rules that rely on group information for each sample.
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
    Rscript {params.rscript1} '{params.outdir}' '{output.outhtml}' \
    '{input.file1}' '{input.file2}' '{params.projectId}' '{params.projDesc}' '{params.rscript2}'
    """
