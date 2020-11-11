# RNA-seek

![Build](https://github.com/skchronicles/RNA-seek/workflows/Tests/badge.svg)  [![GitHub issues](https://img.shields.io/github/issues/skchronicles/RNA-seek)](https://github.com/skchronicles/RNA-seek/issues)  [![GitHub license](https://img.shields.io/github/license/skchronicles/RNA-seek)](https://github.com/skchronicles/RNA-seek/blob/main/LICENSE)

An open-source, reproducible, and scalable solution for analyzing RNA-seq data.

### 1. Getting Started

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (TBA later!).

##### 1.1 Download the workflow
Please [clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/skchronicles/RNA-seek.git
# Change your working directory to the RNA-seek repo
cd RNA-seek/
```
##### 1.2 Add snakemake to PATH
Please make sure that `snakemake>=5.19` is in your `$PATH`. If you are in Biowulf, please load the following environment module:
```bash
# Recommend running snakemake>=5.19
module load snakemake/5.24.1
```

##### 1.3 Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

##### 1.4 Lint the workflow

Please lint the workflow prior to dry-running.
```bash
snakemake --lint -s workflow/Snakefile
```
##### 1.5 Dry-run the workflow

Please run the following command to dry-run the snakemake pipeline:
```bash
snakemake -npr -s workflow/Snakefile
```

### 2. Usage

Submit master job to the cluster:
```bash
# Local Testing (Interactive node)
sinteractive --mem=110g --cpus-per-task=12 --gres=lscratch:200
module load singularity snakemake

## Method 1: Using environment modules
snakemake -pr --use-envmodule --cores 12 --configfile=.tests/run.json

## Method 2: Using Singularity Images
snakemake -npr --use-singularity --singularity-args '-B /data/CCBR_Pipeliner/db/PipeDB/,/lscratch,/fdb' --cores 12 --configfile=.tests/run.json

# Generate a Snakemake Report
module load graphviz # dot needs to be in $PATH
snakemake --report .tests/report.html --cores 12 --configfile=.tests/run.json

# Add later
echo "Coming soon!"
```

### 3. Contribute

This section is for new developers working with the RNA-seek pipeline. If you have added new features or adding new changes, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or org account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local filesystem.
3. Copy the modified files to the cloned fork.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) to this repository.

