## Reference genomes
> _**Warning:**_ This section contains FTP links for downloading each reference file.  

The quality-control and differential expression pipeline support the following genomes:    
<!---
**Human** `hg19` `hg38` `hg38_30` `hs37d5` `hs38d1` `hg38_30`  
**Human + _Integrated Virus_** `hg38_30_KSHV` `hg38_HPV16`  
**Mouse** `mm10` `mm9` `mm10_M21`  
**Dog** `canFam3`  
**Rhesus macaque** `Mmul_8.0.1`
--->

| GenomeName   | Species                                   | Annotation Version                                           | Comments                                                     |
| ------------ | ----------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| hg19         | Homo sapiens (human)                      | [Gencode Release 19](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz) | [GRCh37](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz), Release date: 07/2013 |
| hg38         | Homo sapiens (human)                      | [Gencode Release 28](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz) | [GRCh38](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz), Annotation Release date: 11/2017 |
| hg38_30      | Homo sapiens (human)                      | [Gencode Release 30](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz) | [GRCh38](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz), Annotation Release date: 11/2018 |
| hs37d5       | Homo sapiens (human)                      | [Gencode Release 19](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz) | [hg19 + decoy sequences](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence) |
| hs38d1       | Homo sapiens (human)                      | [Gencode Release 28](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz) | hg38 + [decoy sequences](https://www.ncbi.nlm.nih.gov/assembly/GCA_000786075.2/) |
| hg38_30_KSHV | Homo sapiens + KSHV                       | [Gencode Release 30](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz) (hg38) + [06/2019](https://www.ncbi.nlm.nih.gov/nuccore/NC_009333.1) (KSHV) | hg38 + [NC_009333.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009333.1). Annotation Release dates: 11/2018(human) + 06/2019(KSHV) |
| hg38_HPV16   | Homo sapiens + HPV16                      | [Gencode Release 28](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz) (hg38) + 03/2019 (HPV16 custom annotation from Zheng lab) | hg38 + HPV16 custom sequence based off of KU298885.1 with custom annotation |
| mm9          | Mus musculus (house mouse)                | [M1](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz) | [NCBIM37](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/NCBIM37.genome.fa.gz), Annotation Release date: 12/2011 |
| mm10         | Mus musculus (house mouse)                | [M18](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.annotation.gtf.gz) | [GRCm38](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.primary_assembly.genome.fa.gz), Annotation Release date: 07/2018 |
| mm10_M21     | Mus musculus (house mouse)                | [M21](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz) | [GRCm38](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.primary_assembly.genome.fa.gz), Annotation Release date: 04/2019 |
| canFam3      | Canis lupus familiaris (dog)              | [Ensembl Release 94](ftp://ftp.ensembl.org/pub/release-94/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.94.gtf.gz)                                                            | [CanFam3.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000002285.3/) |
| Mmul_8.0.1   | Macaca mulatta (Rhesus monkey or macaque) | [Ensembl Release 97](ftp://ftp.ensembl.org/pub/release-97/gtf/macaca_mulatta/Macaca_mulatta.Mmul_8.0.1.97.gtf.gz) | [Mmul_8.0.1 (rheMac8)](ftp://ftp.ensembl.org/pub/release-95/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa.gz) |

> **Please note**: If you are looking for a reference genome and/or annotation _**that is currently not available**_, it can be generated using [Pipeliner Index Maker (PIM) ](https://github.com/CCBR/PIM). Given the reference's FASTA file `ref.fa` and a GTF file `genes.gtf`, PIM will create all of the required reference files to run RNA-seq pipeline on Biowulf.

## Tools and versions

### Quality-control pipeline
> _Raw data > Adapter Trimming > Alignment > Quantification (genes and isoforms)_

| Tool           | Version | Notes                                                                                                         |
|----------------|:-------:|---------------------------------------------------------------------------------------------------------------|
| FastQC<sup>2</sup>         |  0.11.5 | **Quality-control step** to assess sequencing quality, run before and after adapter trimming                          |
| Cutadapt<sup>3</sup>       |   1.18  | **Data processing step** to remove adapter sequences and perform quality trimming                             |
| Kraken<sup>18</sup>         |   1.1   | **Quality-control step** to assess microbial taxonomic composition                                            |
| KronaTools<sup>19</sup>     |   2.7   | **Quality-control step** to visualize kraken output                                                           |
| FastQ Screen<sup>21</sup>   |  0.9.3  | **Quality-control step** to assess contamination; additional dependencies: `bowtie2/2.3.4`, `perl/5.24.3`     |
| STAR<sup>4</sup>           |  2.7.0f | **Data processing step** to align reads against reference genome (using its two-pass mode)                    |
| QualiMap<sup>20</sup>       |  2.2.1  | **Quality-control step** to assess various alignment metrics, also calculates insert_size                     |
| Picard<sup>12</sup>         | 2.17.11 | **Quality-control step** to run `MarkDuplicates`, `CollectRnaSeqMetrics` and `AddOrReplaceReadGroups`         |
| Preseq<sup>1</sup>         |  2.0.3  | **Quality-control step** to estimate library complexity                                                       |
| SAMtools<sup>17</sup>       |   1.6   | **Quality-control step** to run `flagstat` to calculate alignment statistics                                  |
| bam2strandedbw | [custom](https://github.com/CCBR/Pipeliner/blob/master/Results-template/Scripts/bam2strandedbw.pe.sh)   | **Summarization step** to convert STAR aligned PE bam file into forward and reverse strand bigwigs suitable for a genomic track viewer like IGV              |
| RSeQC<sup>11</sup>          | 2.6.4   | **Quality-control step** to infer stranded-ness and read distributions over specific genomic features         |
| RSEM<sup>5</sup>           | 1.3.0   | **Data processing step** to quantify gene and isoform counts                                                  |
| Subread<sup>14</sup>        | 1.5.2   | **Data processing step** to run `featureCounts`, an alternative quantification method to RSEM                 |
| PCA Report<sup>16</sup>      | [custom](https://github.com/CCBR/Pipeliner/blob/master/Results-template/Scripts/PcaReport.Rmd)   | **Summarization step** to identify outliers prior to DE, contains pre- and post- normalization plots  |
| MultiQC<sup>15</sup>        | 1.4     | **Reporting step** to aggregate sample statistics and quality-control information across all sample           |

### Differential expression pipeline
> _Raw counts matrix > Normalization > Differential Expression Analysis > Fuctional Impact_

| Tool          |                                              Version                                             | Notes                                                                                                         |
|---------------|:------------------------------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------------------------|
| filtersamples<sup>16</sup> | [custom](https://github.com/CCBR/Pipeliner/blob/master/Results-template/Scripts/filtersamples.R) | **Data processing step** to remove low CPM genes prior to differential expression analysis                    |
| PCAReport<sup>16</sup>     | [custom](https://github.com/CCBR/Pipeliner/blob/master/Results-template/Scripts/PcaReport.Rmd)   | **Summarization step** to identify outliers prior to DE, contains pre- and post- normalization plots  |
| EBSeq<sup>22</sup>         | 1.2.0                                                                                            | **Data processing step** to find differentially expressed isoforms, additional dependencies: `rsem/1.3.0`     |
| edgeR<sup>23</sup>         | 3.24.3                                                                                           | **Data processing step** to find differentially expressed genes. Counts are modeled using a negative binomial distribution with mean equal to the multiplication of library size and relative abundance while a quasi-likelihood F-test is used for testing gene differential expression |
| DESeq2<sup>13</sup>        | 1.22.2                                                                                           | **Data processing step** to find differentially expressed genes. Counts are modeled using a negative binomial distribution similar to edgeR while a wald-test is implemented to test for differential expression |
| limma<sup>7,8</sup>         | 3.38.3                                                                                           | **Data processing step** to find differentially expressed genes. Log-transformed counts are modeled using a method analogous to a t-distribution while a moderated t-statistics is used to test for differential expression                                          |
| l2p<sup>16</sup>     | [custom](https://github.com/CCBR/l2p)   | **Summarization step** for gene set enrichment analysis  |

## References
<sup>**1.**	Daley, T. and A.D. Smith, Predicting the molecular complexity of sequencing libraries. Nat Methods, 2013. 10(4): p. 325-7.</sup>  
<sup>**2.**    Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.</sup>  
<sup>**3.**	Martin, M. (2011). "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet 17(1): 10-12.</sup>  
<sup>**4.**	Dobin, A., et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 2013. 29(1): p. 15-21.</sup>  
<sup>**5.**	Li, B. and C.N. Dewey, RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics, 2011. 12: p. 323.</sup>  
<sup>**6.**	Harrow, J., et al., GENCODE: the reference human genome annotation for The ENCODE Project. Genome Res, 2012. 22(9): p. 1760-74.</sup>  
<sup>**7.**	Law, C.W., et al., voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol, 2014. 15(2): p. R29.</sup>  
<sup>**8.**	Smyth, G.K., Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Stat Appl Genet Mol Biol, 2004. 3: p. Article3.</sup>  
<sup>**9.**	Fabregat, A., et al., The Reactome Pathway Knowledgebase. Nucleic Acids Res, 2018. 46(D1): p. D649-D655.</sup>  
<sup>**10.**	Liberzon, A., et al., Molecular signatures database (MSigDB) 3.0. Bioinformatics, 2011. 27(12): p. 1739-40.</sup>  
<sup>**11.**    Wang, L., et al. (2012). "RSeQC: quality control of RNA-seq experiments." Bioinformatics 28(16): 2184-2185.</sup>  
<sup>**12.**    The Picard toolkit. https://broadinstitute.github.io/picard/.</sup>  
<sup>**13.**    Love, M. I., et al. (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biol 15(12): 550.</sup>  
<sup>**14.**    Liao, Y., et al. (2013). "The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote." Nucleic Acids Research 41(10): e108-e108.</sup>  
<sup>**15.**    Ewels, P., et al. (2016). "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32(19): 3047-3048.</sup>  
<sup>**16.**    R Core Team (2018). R: A Language and Environment for Statistical Computing. Vienna, Austria, R Foundation for Statistical Computing.</sup>  
<sup>**17.**    Li, H., et al. (2009). "The Sequence Alignment/Map format and SAMtools." Bioinformatics 25(16): 2078-2079.</sup>  
<sup>**18.**    Wood, D. E. and S. L. Salzberg (2014). "Kraken: ultrafast metagenomic sequence classification using exact alignments." Genome Biol 15(3): R46.</sup>  
<sup>**19.**    Ondov, B. D., et al. (2011). "Interactive metagenomic visualization in a Web browser." BMC Bioinformatics 12(1): 385.</sup>  
<sup>**20.**    Okonechnikov, K., et al. (2015). "Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data." Bioinformatics 32(2): 292-294.</sup>  
<sup>**21.**    Wingett, S. and S. Andrews (2018). "FastQ Screen: A tool for multi-genome mapping and quality control." F1000Research 7(2): 1338.</sup>  
<sup>**22.**    Leng, N., et al. (2013). "EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments." Bioinformatics 29(8): 1035-1043.</sup>  
<sup>**23.**    Robinson, M. D., et al. (2009). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics 26(1): 139-140.</sup>  
