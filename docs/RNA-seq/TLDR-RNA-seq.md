## 1. Introduction

When processing RNA-sequencing data, there are often many steps that we must repeat. These are usually steps like removing adapter sequences, aligning reads against a reference genome, checking the quality of the data, and quantifying counts. RNA-seek is composed of several sub commands or convience functions to automate these repetitive steps.

With RNA-seek, you can run your samples through our highly-reproducible pipeline, build resources for new reference genomes, and more!

Here is a list of available rna-seek `sub commands`:
   
 - **`run`**: run the rna-seq pipeline   
 - **`build`**: build reference files   
 - **`cache`**: cache remote resources locally  
 - **`unlock`**: unlock a working directory  

> This page contains information for building reference files and running the RNA-seek pipeline. For more information about each of the available sub commands, please see the [usage section](). 

## 2. Setup RNA-seek
_Estimated Reading Time: 3 Mintutes_

RNA-seek has two dependencies: `singularity` and `snakemake`. These dependencies can be installed by a sysadmin; however, snakemake is readily available through conda. Before running the pipeline or any of the commands below, please ensure singularity and snakemake are in your `$PATH`. Please see follow the instructions below for getting started with the RNA-seek pipeline.

### 2.1 Login to cluster
```bash
# Setup Step 0.) ssh into cluster's head node
# example below for Biowulf cluster
ssh -Y $USER@biowulf.nih.gov
```


### 2.2 Grab an interactive node
```bash 
# Setup Step 1.) Please do not run RNA-seek on the head node!
# Grab an interactive node first
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
```

### 2.3 Load dependecies
```bash 
# Setup Step 2.) Add singularity and snakemake executables to $PATH
module purge
module load singularity snakemake

# Setup Step 3.) Download RNA-seek and add to $PATH
# Clone the RNA-seek repository from Github
git clone https://github.com/skchronicles/RNA-seek.git
export PATH=${PWD}/RNA-seek:${PATH}
```

## 3. Building Reference files

In this example, we will start off by building reference files downloaded from [GENCODE](https://www.gencodegenes.org/). We recommend downloading the `PRI` Genome FASTA file and annotation from [GENCODE](https://www.gencodegenes.org/). These `PRI` reference files contain the primary chromosomes and scaffolds. We **do not** recommend downloading the `CHR` reference files! 

Here is more information about GENCODE's [v36 release](https://www.gencodegenes.org/human/release_36.html) for the human reference genome.

### 3.1 Download References from GENCODE

```bash
# Build Step 0.) Please do not run RNA-seek on the head node!
# Grab an interactive node first
# Assumes that you have already ssh-ed into cluster
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash

# Build Step 1.) Download the PRI Genome FASTA file for GRCh38.p13
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
gzip -d GRCh38.primary_assembly.genome.fa.gz

# Build Step 2.) Download the PRI release 36 annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz
gzip -d gencode.v36.primary_assembly.annotation.gtf.gz
```

### 3.2 Run Build pipeline 
```bash
# Build Step 3.) Load dependencies
module purge
module load singularity snakemake

# Build Step 4.) Dry-run the build pipeline
rna-seek build  --ref-fa GRCh38.primary_assembly.genome.fa \
                --ref-name hg38 \
                --ref-gtf gencode.v36.primary_assembly.annotation.gtf \
                --gtf-ver 36 --output /scratch/$USER/hg38_36 --dry-run


# Build Step 5.) Submit the build pipeline to cluster
rna-seek build  --ref-fa GRCh38.primary_assembly.genome.fa \
                --ref-name hg38 \
                --ref-gtf gencode.v36.primary_assembly.annotation.gtf \
                --gtf-ver 36 --output /scratch/$USER/hg38_36 
```

An email notification will be sent out when the pipeline starts and ends. Once the build pipeline completes, you can run RNA-seek with the provided test dataset. Please see the intructions below for more information.

## 4. Running RNA-seek 

Run RNA-seek with the reference files we built above using hg38 (GRCh38.p13) Genome FASTA file and GENCODE release 36 annotation (GTF). For more information about how the reference files we generated, please see the intructions above. You can use those instructions as a guide for building any new reference genomes in the future. 


### 4.1 Dry-run pipeline 

Dry-run the pipeline prior to submiting the pipeline's master job. Please note that if you wish to run RNA-seek with a new dataset, you will only need to update the values provided to the `--input` and `--output` arguments (and maybe `--genome`). The `--input` argument supports globbing. If this is the first time running RNA-seek with for given dataset, the `--output` directory should _**not**_ exist on your local filesystem. It will be created automatically during runtime.

```bash
# Run Step 0.) Please do not run RNA-seek on the head node!
# Grab an interactive node first
# Assumes that you have already ssh-ed into cluster
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash

# Run Step 1.) Load dependencies
module purge
module load singularity snakemake

# Run Step 2.) Dry-run the pipeline with test dataset
# And reference genome generated in the steps above
# Test data consists of sub sampled FastQ files 
rna-seek run \
    --input RNA-seek/.tests/*.R?.fastq.gz \
    --output /scratch/${USER}/runner_hg38_36/  \
    --genome /scratch/${USER}hg38_36/hg38_36.json  \
    --mode slurm \
    --star-2-pass-basic \
    --dry-run
```

### 4.2 Run pipeline 

Kick off the pipeline by submiting the master job to the cluster. It is essentially the same command above without the `--dry-run` flag. 

```bash
# Run Step 3.) Submit the master job
# Runs the RNA-seek pipeline  with the 
# reference genome generated in the steps above
# and with the test dataset
rna-seek run \
    --input RNA-seek/.tests/*.R?.fastq.gz \
    --output /scratch/${USER}/runner_hg38_36/  \
    --genome /scratch/${USER}/hg38_36/hg38_36.json  \
    --mode slurm \
    --star-2-pass-basic \
    --dry-run
```

An email notification will be sent out when the pipeline starts and ends.