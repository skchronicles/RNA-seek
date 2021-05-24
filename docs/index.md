# RNA-seek

![Docker Pulls](https://img.shields.io/docker/pulls/nciccbr/ccbr_arriba_2.0.0) [![Build](https://github.com/skchronicles/RNA-seek/workflows/Tests/badge.svg)](https://github.com/skchronicles/RNA-seek/actions)  [![GitHub issues](https://img.shields.io/github/issues/skchronicles/RNA-seek?color=brightgreen)](https://github.com/skchronicles/RNA-seek/issues)  [![GitHub license](https://img.shields.io/github/license/skchronicles/RNA-seek)](https://github.com/skchronicles/RNA-seek/blob/main/LICENSE) 

An open-source, reproducible, and scalable solution for analyzing RNA-sequencing data.

## 1. Introduction
Welcome to RNA-seek's documentation! 

This guide is the main source of documentation for users that are getting started with the [RNA-seek pipeline](https://github.com/skchronicles/RNA-seek). If you are not familiar with RNA-sequencing, please checkout our [theory and practical guide](RNA-seq/Theory.md). That section provides a conceptual overview to RNA-seq analysis and as well as a set of generalized guidelines to interpret different quality-control metrics.  If you are a new user, we highly recommend reading through our [getting started](RNA-seq/TLDR-RNA-seq.md) section. This page contains information needed to quickly build new reference files and setup the pipeline for running in your compute environment. 

RNA-seek is composed several inter-related sub commands to faciliate the analysis of RNA-sequencing data. For more information about each available sub command, please see the usage section. To help out new users, an example of each command is also provided. The [resources page](RNA-seq/Differential-expression-pipeline-tools-and-versions.md) contains more information about the pipeline's default reference genomes along with every tool and Docker image the pipeline employs. 

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](troubleshooting.md) prior to [opening an issue on Github](https://github.com/skchronicles/RNA-seek/issues).

## 2. Overview

**RNA-seek** is a comprehensive, open-source RNA-seq pipeline that relies on technologies like [Docker<sup>1</sup>](https://www.docker.com/why-docker) and [Singularity<sup>2</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>3</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster or cloud provider.

RNA-seek can be run locally on a compute instance, on-premise using a cluster, or on the cloud using AWS. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM, or run on AWS using Tibanna (feature coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

A bioinformatics pipeline is more than the sum of its data processing steps. A pipeline without quality-control steps provides a myopic view of the potential sources of variation within your data (i.e., biological verses technical sources of variation). RNA-seek pipeline is composed of a series of quality-control and data processing steps.

![RNA-seq quantification pipeline](RNA-seq/images/RNA-seek_Pipeline.svg){width="650"; align=right} <br><sup>**Fig 1. An Overview of the Quantification and Quality-control Pipeline.** Gene and isoform counts are quantified and a series of QC-checks are performed to assess the quality of the data. This pipeline stops at the generation of a raw counts matrix, which is input to the next sub-workflow. To run the pipeline, a user must select their raw data directory (i.e. the location to their FastQ files), a reference genome, and output directory (i.e. the location where the pipeline performs the analysis). Quality-control information is summarized across all samples in the MultiQC and RNA report.</sup>

## 3. Pipeline

The accuracy of the downstream interpretations made from transcriptomic data are highly dependent on initial sample library. Unwanted sources of technical variation, which if not accounted for properly, can influence the results. 

In addition to generating a MultiQC report, the RNA-seek pipeline also generates a [rNA Report](https://github.com/CCBR/rNA) to interactively allow users to identify problematic samples prior to performing any downstream analysis. RNA-seek's comprehensive quality-control helps ensure your results are reliable and reproducible across experiments.  In the data processing steps, RNA-seek quantifies gene and isoform expression and predicts gene fusions. Please note that the detection of alternative splicing events and variant calling will be incorporated in a later release.


### 3.1 Quality Control
**Quality Control**   
[*FastQC*<sup>4</sup>](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used to assess the sequencing quality. FastQC is run twice, before and after adapter trimming. It generates a set of basic statistics to identify problems that can arise during sequencing or library preparation. FastQC will summarize per base and per read QC metrics such as quality scores and GC content. It will also summarize the distribution of sequence lengths and will report the presence of adapter sequences.
 
[*Kraken2*<sup>5</sup>](http://ccb.jhu.edu/software/kraken2/) and [*FastQ Screen*<sup>6</sup>](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) are used to screen for various sources of contamination. During the process of sample collection to library preparation, there is a risk for introducing wanted sources of DNA. FastQ Screen compares your sequencing data to a set of different reference genomes to determine if there is contamination. It allows a user to see if the composition of your library matches what you expect. Also, if there are high levels of microbial contamination, Kraken can provide an estimation of the taxonomic composition. Kraken can be used in conjunction with [*Krona*<sup>7</sup>](https://github.com/marbl/Krona/wiki/KronaTools) to produce interactive reports.

[*Preseq*<sup>8</sup>](http://smithlabresearch.org/software/preseq/) is used to estimate the complexity of a library for each samples. If the duplication rate is very high, the overall library complexity will be low. Low library complexity could signal an issue with library preparation where very little input RNA was over-amplified or the sample may be degraded.

[*Picard*<sup>9</sup>](https://broadinstitute.github.io/picard/) can be used to estimate the duplication rate, and it has another particularly useful sub-command called CollectRNAseqMetrics which reports the number and percentage of reads that align to various regions: such as coding, intronic, UTR, intergenic and ribosomal regions. This is particularly useful as you would expect a library constructed with ploy(A)-selection to have a high percentage of reads that map to coding regions. Picard CollectRNAseqMetrics will also report the uniformity of coverage across all genes, which is useful for determining whether a sample has a 3' bias (observed in ploy(A)-selection libraries containing degraded RNA).

[*RSeQC*<sup>10</sup>](http://rseqc.sourceforge.net/) is another particularity useful package that is tailored for RNA-seq data. It is used to calculate the inner distance between paired-end reads and calculate TIN values for a set of canonical protein-coding transcripts. A median TIN value is calucated for each sample, which analogous to a computationally derived RIN.

[MultiQC<sup>11</sup>](https://multiqc.info/) is used to aggreate the results of each tool into a single interactive report.  

### 3.2 Data Processing
 
[*Cutadapt*<sup>12</sup>](https://cutadapt.readthedocs.io/en/stable/) is used to remove adapter sequences, perform quality trimming, and remove very short sequences that would otherwise multi-map all over the genome prior to alignment. 

[*STAR*<sup>13</sup>](https://github.com/alexdobin/STAR) is used to align reads to the reference genome. The RNA-seek pipeline runs STAR in a two-passes where splice-junctions are collected and aggregated across all samples and provided to the second-pass of STAR. In the second pass of STAR, the splice-junctions detected in the first pass are inserted into the genome indices prior to alignment.

[*RSEM*<sup>14</sup>](https://github.com/deweylab/RSEM) is used to quantify gene and isoform expression. The expected counts from RSEM are merged across samples to create a two counts matrices for gene counts and isoform counts.

[*Arriba*<sup>15</sup>](https://arriba.readthedocs.io/en/latest/) is used to predict gene-fusion events. The pre-built human and mouse reference genomes use Arriba blacklists to reduce the false-positive rate.

## 4. Contribute

This site is a living document, created for and by the genomics community. RNA-seek is maintained by the developers and bioinformaticians at the NIH and is improved by feedback from external collaborators like *you*! 

We want to make it easy for users to connect with us to share ideas, solve problems, and to continuously deliver the best pipelines. We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/CCBR/RNA-seek). Please checkout our [project's roadmap](https://github.com/skchronicles/RNA-seek/projects/1) to see our planned improvements for the next release of RNA-seek.

## 5. References
<sup>**1.**  Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2.</sup>  
<sup>**2.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**3.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup> 
<sup>**4.** Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.</sup>  
<sup>**5.** Wood, D. E. and S. L. Salzberg (2014). "Kraken: ultrafast metagenomic sequence classification using exact alignments." Genome Biol 15(3): R46.</sup>  
<sup>**6.** Wingett, S. and S. Andrews (2018). "FastQ Screen: A tool for multi-genome mapping and quality control." F1000Research 7(2): 1338.</sup>  
<sup>**7.** Ondov, B. D., et al. (2011). "Interactive metagenomic visualization in a Web browser." BMC Bioinformatics 12(1): 385.</sup>  
<sup>**8.** Daley, T. and A.D. Smith, Predicting the molecular complexity of sequencing libraries. Nat Methods, 2013. 10(4): p. 325-7.</sup>  
<sup>**9.** The Picard toolkit. https://broadinstitute.github.io/picard/.</sup>  
<sup>**10.** Wang, L., et al. (2012). "RSeQC: quality control of RNA-seq experiments." Bioinformatics 28(16): 2184-2185.</sup>  
<sup>**11.** Ewels, P., et al. (2016). "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32(19): 3047-3048.</sup>  
<sup>**12.** Martin, M. (2011). "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet 17(1): 10-12.</sup>  
<sup>**13.** Dobin, A., et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 2013. 29(1): p. 15-21.</sup>  
<sup>**14.** Li, B. and C.N. Dewey, RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics, 2011. 12: p. 323.</sup>  
<sup>**15.** Uhrig, S., et al. (2021). "Accurate and efficient detection of gene fusions from RNA sequencing data". Genome Res. 31(3): 448-460.</sup>  




<!-- Relative links -->
  [1]: contact-us.md
