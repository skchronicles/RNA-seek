import json
input_fa="/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/mm10.fa"
input_gtf="/data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/mm10/gencode.vM21.annotation.gtf"
params_workdir="/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_M21"
params_genome="mm10_M21"
output_json=params_workdir+"/"+params_genome+".json"
bigdict=dict()
bigdict["references"]=dict()
for i in ["exomeseq", "genomeseq", "rnaseq", "rnaseqvargerm", "ChIPseq"]:
    bigdict["references"][i]=dict()
bigdict["references"]["rnaseq"]["GENOMEFILE"]=input_fa
bigdict["references"]["rnaseq"]["GENOME"]=input_fa
bigdict["references"]["rnaseq"]["GTFFILE"]=input_gtf
bigdict["references"]["rnaseq"]["STARDIR"]=params_workdir+"/STAR/2.7.0f/genes-"
bigdict["references"]["rnaseq"]["STARREF"]=params_workdir+"/STAR/2.7.0f/genes-"
bigdict["references"]["rnaseq"]["ANNOTATE"]=params_workdir+"/annotate.genes.txt"
bigdict["references"]["rnaseq"]["ANNOTATEISOFORMS"]=params_workdir+"/annotate.isoforms.txt"
bigdict["references"]["rnaseq"]["REFFLAT"]=params_workdir+"/refFlat.txt"
bigdict["references"]["rnaseq"]["BEDREF"]=params_workdir+"/genes.ref.bed"
bigdict["references"]["rnaseq"]["GENEINFO"]=params_workdir+"/geneinfo.bed"
bigdict["references"]["rnaseq"]["KARYOBEDS"]=params_workdir+"/karyobeds"
bigdict["references"]["rnaseq"]["RSEMREF"]=params_workdir+"/rsemref/"+params_genome
bigdict["references"]["rnaseq"]["RRNALIST"]=params_workdir+"/"+params_genome+".rRNA_interval_list"
bigdict["references"]["rnaseq"]["FASTQ_SCREEN_CONFIG"]="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf"
bigdict["references"]["rnaseq"]["FASTAWITHADAPTERSETC"]="/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters_new.fa"
bigdict["references"]["rnaseq"]["adapter.file"]="/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters.ngsqc.dat"
bigdict["references"]["rnaseq"]["trimmomatic.adapters"]="/data/CCBR_Pipeliner/db/PipeDB/dev/adapters2.fa"
bigdict["references"]["rnaseq"]["fastqc.adapters"]="/data/CCBR_Pipeliner/db/PipeDB/dev/fastqc.adapters"
bigdict["references"]["rnaseq"]["ORGANISM"]="MOUSE"
with open(output_json, 'w') as fp:
    json.dump(bigdict, fp, indent=4)
