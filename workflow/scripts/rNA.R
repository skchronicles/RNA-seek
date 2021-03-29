#!/usr/bin/env Rscript

# USAGE: Rscript rNA.R -m src/rNA.Rmd -r data/TCGA-GBM_Raw_RSEM_Genes.txt -t data/TCGA-GBM_TINs.txt -q data/multiqc_matrix.txt -o "$PWD"

library(argparse)

# Parse Args
parser <- ArgumentParser()

# rNA Rmarkdown file
parser$add_argument("-m", "--rmarkdown", type="character", required=TRUE,
                    help="Required Input File: rNA Rmarkdown file")

# Raw Counts Matrix
parser$add_argument("-r", "--raw_counts", type="character", required=TRUE,
                    help="Required Input File: Raw counts matrix")

# TIN Counts Matrix
parser$add_argument("-t", "--tin_counts", type="character", required=TRUE,
                    help="Required Input File: TIN counts matrix")

# QC Metadata Table
parser$add_argument("-q", "--qc_table", type="character", required=TRUE,
                    help="Required Input File: QC Metadata Table")

# Output directory
parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
                    help="Required Path to Output Directory")

# Output HTML Filename
parser$add_argument("-f", "--output_filename", type="character", required=FALSE, default = 'rNA.html',
                    help="Optional Output HTML Filename: Defaults to 'rNA.html'")

# Display sample names
parser$add_argument("-a", "--annotate", action="store_true", default=FALSE,
    help="Display sample names in complex heatmap: Defaults to FALSE")


args <- parser$parse_args()

# Get current working directory to setwd() of rmarkdown, else PATHs must be absolute
working_directory = getwd()  # allows for relative paths to this main entry script

# Generate HTML output
rmarkdown::render(args$rmarkdown, output_file=file.path(args$output_dir, args$output_filename), params = list(
  raw = args$raw_counts,
  tin = args$tin_counts,
  qc = args$qc_table,
  wdir = working_directory,
  annot = args$annotate
 )
)
