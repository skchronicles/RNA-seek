# Example Usage: Rscript pcacall.R 'DEG_cntrl-test_0.5_2' 'outfilename.html'  'STAR_files/sampletable.txt' 'DEG_cntrl-test_0.5_2/RawCountFile_RSEM_genes_filtered.txt' 'hg19seDEG' 'Enter CCBR Project Description and Notes here.' '/path/to/workingDir/Scripts/PcaReport.Rmd'
## grab args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
outHtml <- args[2]
pcaRmd <- args[7]
Sys.setenv(RSTUDIO_PANDOC="/usr/local/apps/rstudio/rstudio-1.1.447/bin/pandoc/")
setwd(DIR) # new 
rmarkdown::render(pcaRmd,output_file=outHtml, params = list(
    folder = args[1],
    sampleinfo = args[3],
    data = args[4],
    projectId = args[5],
    projectDesc = args[6]
  )
)
