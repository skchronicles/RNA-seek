## 1. Introduction
RNA-sequencing (*RNA-seq*) has a wide variety of applications; this transcriptome profiling method can be used to quantify gene and isoform expression, find changes in alternative splicing, detect gene-fusion events, call variants and much more. 

It is also worth noting that RNA-seq can be coupled with other biochemical assays to analyze many other aspects of RNA biology, such as RNA–protein binding (CLIP-seq, RIP-seq), RNA structure (SHAPE-seq), or RNA–RNA interactions (CLASH-seq). These applications are, however, beyond the scope of this documentation as we focus on *typical* RNA-seq project (i.e. quantifying expression and gene fusions). Our focus is to outline current standards and resources for the bioinformatics analysis of RNA-seq data. We do not aim to provide an exhaustive compilation of resources or software tools. Rather, we aim to provide a guideline and conceptual overview for RNA-seq data analysis based on our best-practices RNA-seq pipeline.

Here we review all of the *typical* major steps in RNA-seq data analysis, starting from experimental design, quality control, read alignment, quantification of gene and transcript levels, and visualization.

## 2. Experimental Design
Just like any other scientific experiment, a good RNA-seq experiment is hypothesis-driven. If you cannot describe the problem you are trying to address, throwing NGS at the problem is not a cure-all solution. Fishing for results is a waste of your time and is bad science. As so, designing a well-thought-out experiment around a testable question will maximize the likelihood of generating high-impact results.

The data that is generated will determine whether you have the potential to answer your biological question of interest. As a prerequisite, you need to think about how you will construct your libraries; the correct sequencing depth to address your question of interest; the number of replicates, and strategies to reduce/mitigate batch effects.

### 2.1 Library construction
rRNA can comprise up to 80% of the RNA in a cell. An important consideration is the RNA extraction protocol that will be used to remove the highly abundant ribosomal RNA (rRNA). For eukaryotic cells, there are two major considerations: choosing whether to enrich for mRNA or whether to deplete rRNA. 

#### 2.1.1 mRNA
Poly-(A) selection is a common method used to enrich for mRNA. This method generates the highest percentage of reads which will ultimately map to protein-coding genes-- making it a common choice for most applications. That being said, poly(A)-selection requires your RNA to be of high quality with minimal degradation. Degraded samples that are followed with ploy(A)-selection may result in a 3’ bias, which in effect, may introduce downstream biases into your results. 

#### 2.1.2 total RNA
The second method captures total RNA through the depletion of rRNA. This method allows you to examine both mRNA and other non-coding RNA species such as lncRNAs. Again, depending on the question you are trying to answer this may be the right method for you. Although, it should be noted that both methods, mRNA and total RNA, require RINs (>8). But if you samples do contain slightly degraded RNA, you might be able to use the total RNA method over poly(A)-selection. 

### 2.2 Sequencing Depth
Sequencing depth or library size is another important design factor. As sequencing depth is increased, more transcripts will be detected (up until a saturation point), and their relative abundance will be quantified more accurately.

At the end of the day, the targeted sequencing depth depends on the aims of the experiment. Are you trying to quantify differences in gene expression, are you trying to quantify differential isoform usage or alternative splicing events? The numbers quoted below are more or less tailored to quantify differences in gene expression. If you are trying to quantify changes in alternative splicing or isoform regulation, you are going to much higher coverage (~ 100M paired-end reads). 

#### 2.2.1 mRNA
For mRNA libraries or libraries generated from a prep kit using poly-(A) selection, we recommend a minimum sequencing depth of 10-20M paired-end reads (or 20-40M reads). RNA must be of high quality or a 3' bias may be observed. 


#### 2.2.2 total RNA
For total RNA libraries, we recommend a sequencing depth of 25-60M paired-end reads (or 50-120M reads). RNA must be of high quality.

> Note: In the sections above and below, when I say to paired-end reads I am referring to read pairs generated from paired-end sequencing of a given cDNA fragment. You will sometimes see reads reported as pairs of reads or total reads.

### 2.3 Replicates
Another important design factor is the number of replicates. That being said, biological replicates are always preferred over technical replicates. 

#### 2.3.1 Recommendation
We recommend 4 biological replicates per experimental condition or group. Having more replicates is good for several reasons because in the real world problems arise. If you have a bad sample that cannot be used due to severe QC issues, you are still left with 3 biological replicates.  This allows you to drop a bad sample without comprising statistical power downstream. 

#### 2.3.2 Bare Minimum
If cost is a factor, at a minimum, 3 biological replicates will ensure good statistical power for downstream analysis. 

### 2.4 Reducing Batch Effects 

Batch effects represent unwanted sources of technical variation. Batch effects introduce non-biological variation into your data, which if not accounted for can influence the results. Through the process of library preparation to sequencing, there are a number of steps (such as RNA extraction to adapter ligation to lane loading, etc.) that might introduce biases into the resulting data. 

As a general rule of thumb, the best way to reduce the introduction of batch effects is through uniform processing-- meaning you need to ensure that differences in sample handling are minimal. This means that samples should be processed by the same lab technician and everything should be done in a uniform manner. That being said, do not extract your RNA at different times, do not use different lots of reagents! If a large number of samples are being processed and everything cannot be done at the same time, process representative samples from each biological group at the same time. This will ensure that batches and your variable of interest do not become confounded. Also, keep note of which samples belong to each batch. This information will be needed for batch correction. 

To reduce the possibility of introducing batch effects from sequencing, all samples should be multiplexed together on the same lane(s). 

| Sample   | Group | **Batch** | <b>Batch*</b> |
|----------|:-----:|:-----:|:------:|
| Treatment_rep_1    |   KO  |   1   |    1   |
| Treatment_rep_2    |   KO  |   2   |    1   |
| Treatment_rep_3    |   KO  |   1   |    1   |
| Treatment_rep_4    |   KO  |   2   |    1   |
| Control_rep_1 |   WT  |   1   |    2   |
| Control_rep_2 |   WT  |   2   |    2   |
| Control_rep_3 |   WT  |   1   |    2   |
| Control_rep_4 |   WT  |   2   |    2   |

> **Batch** = properly balanced batches, easily corrected :relaxed:  
> <b>Batch*</b> = groups and batch totally confounded, cannot be corrected :worried:

That being said, some problems cannot be bioinformatically corrected. If your variable of interest is totally confounded with your batches, applying batch correction to fix the problem is not going to work, and will lead to undesired results (i.e. `Batch*` column). If batches must be introduced due to other constraining factors, please keep note which samples belong to each batch, and please put some thought into how to properly balance samples across your batches.

## 3. Quality Control
Quality-control (**QC**) is extremely important! As the old adage goes: *Garbage in, Garbage out!* If there is one thing that to take away from this document, let it be that. Performing QC checks will help ensure that your results are reliable and reproducible.

It is worth noting that there is a large variety of open-source tools that can be used to assess the quality of your data so there is no reason to re-invent the wheel. Please keep this in mind but also be aware that there are many wheels *per se*, and you will need to know which to use and when. In this next section, we will cover different quality-control checks that can be applied at different stages of your RNA-seq analysis. These recommendations are based on a few tools our best-practices RNA-seq pipeline employs.

### 3.1 Pre-aligment
Before drawing biological conclusions, it is important to perform quality control checks to ensure that there are no signs of sequencing error, biases in your data, or other sources of contamination. Modern high-throughput sequencers generate millions of reads per run, and in the real world, problems can arise. 

The general idea is to assess the quality of your reads before and after adapter removal and to check for different sources of contamination before proceeding to alignment. Here are a few of the tools that we use and recommend.

#### 3.1.1 Sequencing Quality
To assess the sequencing quality of your data, we recommend running FastQC before and after adapter trimming. FastQC generates a set of basic statistics to identify problems that can arise during sequencing or library preparation. FastQC will summarize per base and per read QC metrics such as quality scores and GC content (ideally, this plot should have a normal distribution with no forms of bimodality). It will also summarize the distribution of sequence lengths and will report the presence of adapter sequences, which is one reason we run it after removing adapters.

#### 3.1.2 Contamination Screening
During the process of sample collection to library preparation, there is a risk for introducing wanted sources of DNA. FastQ Screen compares your sequencing data to a set of different reference genomes to determine if there is contamination. It allows a user to see if the composition of your library matches what you expect. If your data has high levels of human, mouse, fungi, or bacterial contamination, FastQ Screen will tell you. FastQ Screen will tell you what percentage of your library aligns against different reference genomes.

If there are high levels of microbial contamination, Kraken will provide an estimation of the taxonomic composition. Kraken can be used in conjunction with Krona to produce interactive reports.

> Note: Due to high levels of homology between organisms, there may be a small portion of your reads that align to an unexpected reference genome. Again, this should be a minimal percentage of your reads.

### 3.2 Post-alignment
Again, there are many tools available to assess the quality of your data post-alignment, and as stated before, there is no need to re-invent the wheel. Please see the table below for a generalized set of guidelines for different pre/post QC metrics.

#### 3.2.1 Library Complexity 
Preseq can be used to estimate the complexity of a library for each of your samples. If the duplication rate is very high, the overall library complexity will be low. Low library complexity could signal an issue with library preparation or sample preparation (FFPE samples) where very little input RNA was over-amplified or the sample may be degraded.

#### 3.2.2 Library Composition 
Picard has a particularly useful sub-command called CollectRNAseqMetrics which reports the number and percentage of reads that align to various regions: such as coding, intronic, UTR, intergenic and ribosomal regions. This is particularly useful as you would expect a library constructed with ploy(A)-selection to have a high percentage of reads that map to coding regions. Picard CollectRNAseqMetrics will also report the uniformity of coverage across all genes, which is useful for determining whether a sample has a 3' bias (observed in libraries containing degraded RNA).

#### 3.2.3 RNA Quality
This is another particularity useful package that is tailored for RNA-seq data. The package is made up of over 20 sub-module that can be used to do things like calculate the average insert size between paired-end reads (which is useful for GEO upload), annotate the percentage of reads spanning known or novel splice junctions, convert a BAM file into a normalized BigWig file, and infer RNA quality.

### 3.3 Guidelines  
Here is a set of generalized guidelines for different QC metrics. Some of these metrics will vary genome-to-genome depending on the quality of the assembly and annotation but that has been taken into consideration for our set of supported reference genomes. 

| QC Metric Guidelines              |             **mRNA**       |    **total RNA**         |
|-----------------------------------|:--------------------------:|:------------------------:|
| *RNA Type(s)*                     |            Coding          |     Coding + non-coding  |
| *RIN*                             |  >= 8 [low RIN ~ 3' bias]  |             >= 8         |
| *Single-end vs Paired-end*        |          Paired-end        |          Paired-end      |
| *Sequencing Depth*                |        10-20M PE reads     |        25-60M PE reads   |
| *FastQC*                          |           Q30 > 70%        |           Q30 > 70%      |
| *Percent Aligned to Reference*    |             > 70%          |             > 65%        |
| *Million Reads Aligned Reference* |         > 7M PE reads      |       > 16.5M PE reads   |
| *Percent Aligned to rRNA*         |             < 5%           |             < 15%        |
| *Picard RNAseqMetrics*            |         Coding > 50%       |         Coding > 35%     |
| *Picard RNAseqMetrics*            | Intronic + Intergenic < 25%  | Intronic + Intergenic < 40%  |

## 4. Data Processing

Starting from raw data (FastQ files), how do we get a raw counts matrix, or how do we get a list of differential expressed genes? Before feeding your data into an R package for differential expression analysis, it needs to be processed to add biological context to it. In this section, we will talk about the data processing pipeline in more detail-- more specifically focusing on primary and secondary analysis.   

### 4.1 Primary Analysis

>Raw data > Adapter Trimming > Alignment > Quantification

#### 4.1.1 Adapter Trimming  
One of the first steps in this process is to remove any unwanted adapters sequences from your reads in before alignment. Adapters are composed of synthetic sequences and should be removed prior to alignment. Adapter removal is especially important in certain protocols, such as miRNA-seq. When smaller fragments are sequenced it is almost certain there will be some form of adapter contamination.

#### 4.1.2 Alignment  
In the alignment step, we add biological context to the raw data. In this step, we align reads to the reference genome to find where the sequenced fragments originate. 

Accurate alignment of the cDNA fragments (which are derived from RNA) is difficult. Alternative splicing introduces the problem of aligning to non-contiguous regions, and using traditional genomic alignment algorithms can produce inaccurate or low-quality alignments due to the combination of alternative splicing and genomic variation (substitutions, insertions, and deletions). This has lead to the development of *splice-aware* aligners like STAR, which are designed to overcome these issues. STAR can also be run in a two-pass mode for enhanced detection of reads mapping to novel splice junctions. 

#### 4.1.3 Quantification
In the quantification step, the number of reads that mapped to a particular genomic feature (such as a gene or isoform) is counted. It is important to keep in mind that raw counts are biased by a number of factors such as library size, feature-length, and other compositional biases. As so, it is important to normalize your data to remove these biases before summarizing differences between groups of samples. 
