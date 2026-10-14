<div align="center">

  <h1>RNA-seek 🔬</h1>

  **_A comprehensive bulk RNA-seq pipeline_**

  [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4901427-blue)](https://doi.org/10.5281/zenodo.4901427) [![GitHub release](https://img.shields.io/github/v/release/skchronicles/RNA-seek?color=blue&include_prereleases)](https://github.com/skchronicles/RNA-seek/releases) ![Docker Pulls](https://img.shields.io/docker/pulls/nciccbr/ccbr_arriba_2.0.0)<br>
  [![Tests](https://github.com/skchronicles/RNA-seek/actions/workflows/main.yaml/badge.svg)](https://github.com/skchronicles/RNA-seek/actions/workflows/main.yaml) [![Docs](https://github.com/skchronicles/RNA-seek/actions/workflows/docs.yml/badge.svg)](https://github.com/skchronicles/RNA-seek/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/skchronicles/RNA-seek?color=brightgreen)](https://github.com/skchronicles/RNA-seek/issues) [![GitHub license](https://img.shields.io/github/license/skchronicles/RNA-seek)](https://github.com/skchronicles/RNA-seek/blob/main/LICENSE)

  <i>
    Its long-term goals: to assess data quality, quantify gene/isoform expression, and perform differential analyses.
  </i>

</div>

## Overview

Welcome to RNA-seek! Before getting started, we recommend reading through [RNA-seek's documentation](https://openomics.github.io/RNA-seek/).

The **`./rna-seek`** pipeline provides several related sub commands to set up and run the workflow across different systems. Each of the available sub commands perform different functions:

- [<code>rna-seek <b>run</b></code>](https://openomics.github.io/RNA-seek/RNA-seq/run/): Runs the RNA-seek pipeline with your input FastQ files.
- [<code>rna-seek <b>unlock</b></code>](https://openomics.github.io/RNA-seek/RNA-seq/unlock/): Unlocks a previous runs output directory.
- [<code>rna-seek <b>build</b></code>](https://openomics.github.io/RNA-seek/RNA-seq/build/): Builds reference files from a genomic FASTA file and a GTF annotation.
- [<code>rna-seek <b>cache</b></code>](https://openomics.github.io/RNA-seek/RNA-seq/cache/): Caches software containers locally for reuse across runs.

RNA-seek is a comprehensive pipeline for bulk RNA-seq analysis that combines raw data processing, quality control, alignment, and quantification of gene and isoform expression. It also optionally supports gene-fusion detection and differential expression analysis. [Singularity<sup>1</sup>](https://docs.sylabs.io/) runs versioned software containers, while [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/) orchestrates the workflow. The pipeline integrates quality checks throughout the workflow to assess sequencing quality, potential contamination, library complexity, and transcript integrity/degradation. MultiQC summarizes results across samples, helping users identify technical sources of variation and samples that require closer review before proceeding with downstream analyses.

As input, RNA-seek accepts gzipped single-end or paired-end FASTQ files, a reference genome, and an output directory. Single-end and paired-end datasets should be processed in separate runs/output directories. The pipeline can run locally on a compute node with `--mode local` or submit jobs to a SLURM cluster with `--mode slurm` (the default). On Biowulf/Skyline, there are pre-built human and mouse reference genomes; however, new reference genomes can be built on other systems via the `rna-seek build` sub command. Please see the [run](https://openomics.github.io/RNA-seek/RNA-seq/run/) and [build](https://openomics.github.io/RNA-seek/RNA-seq/build/) for more information and examples. 

For help with an issue, check the [FAQ](https://openomics.github.io/RNA-seek/general-questions/) and [troubleshooting guide](https://openomics.github.io/RNA-seek/troubleshooting/) before [opening an issue on GitHub](https://github.com/skchronicles/RNA-seek/issues).

## Dependencies

**Requires:** `singularity>=3.5` and `snakemake<=7.32.4`.

[Snakemake](https://snakemake.readthedocs.io/en/v7.32.4/getting_started/installation.html) and [Singularity](https://docs.sylabs.io/) must be installed on the target system and available in your `PATH`. The command-line entry point also requires python3. To ensure full compatibility, please use a Snakemake version less than or equal to 7.32.4.

Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline relies on versioned docker images from [DockerHub](https://hub.docker.com/). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution.

## Installation

Clone the repository to your local filesystem and inspect the available commands:

```bash
# Clone the repository from GitHub
git clone https://github.com/skchronicles/RNA-seek.git

# Change your working directory
cd RNA-seek/

# Add dependencies to PATH, on Biowulf,
# load these available environment modules
module load snakemake singularity

# Get usage information
./rna-seek -h       # view all available sub commands
./rna-seek run -h   # view options for the running the pipeline
```

On Biowulf, RNA-seek is also available as a system-wide module. Run `module avail openomics/rna-seek` to list available versions of the pipeline. On other systems, load or install snakemake and singularity using your environment’s supported method before running the pipeline.

Follow the [getting started guide](https://openomics.github.io/RNA-seek/RNA-seq/TLDR-RNA-seq/) to prepare references and run the pipeline.

## Contribute

RNA-seek is part of the [OpenOmics](https://openomics.github.io/) pipeline collection and improves through community feedback. We welcome bug reports, documentation improvements, and new features. To contribute, fork the repository, make your changes on a branch, and open a pull request to the [GitHub repository](https://github.com/skchronicles/RNA-seek).

## Cite

If you use RNA-seek in your work, cite the release you used. Release-specific citations are available from the project's [Zenodo](https://zenodo.org/records/4901427). A general software citation is provided below.

<details>
  <summary><b><i>BibTeX</i></b></summary>

```bibtex
@software{openomics_2024_4901427,
  author       = {Kuhn, Skyler and
                  Tandon, Mayank and
                  Routsong, Ryan},
  title        = {skchronicles/RNA-seek},
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.4901427},
  url          = {https://doi.org/10.5281/zenodo.4901427},
}
```

</details>

<details>
  <summary><b><i>APA</i></b></summary>

```text
Kuhn, S., Tandon, M., & Routsong, R. (2024). skchronicles/RNA-seek [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.4901427
```

</details>

Please also cite the tools used in your analysis. See the [resources page](https://openomics.github.io/RNA-seek/RNA-seq/Resources/) for tool descriptions/versions and more information about the default pre-built reference files.

## References

<sup>**1.** Kurtzer, G. M., Sochat, V., and Bauer, M. W. (2017). Singularity: Scientific containers for mobility of compute. PLOS ONE, 12(5), e0177459. [doi:10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459).</sup>

<sup>**2.** Köster, J., and Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520–2522. [doi:10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480).</sup>
