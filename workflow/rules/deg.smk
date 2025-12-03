# Differential gene expression snakemake rules imported in the main Snakefile.
import textwrap
from scripts.common import (
    allocated
)


# Differential expression rules
rule limma_diff_gene_expression:
    """
    Data-processing step to find differentially expressed genes using the R package,
    limma. This step uses the groups and contrasts file provided by the user. The
    groups file is a tab-delimited file containing at least two of the following
    columns:
        1.  Sample: Base name of each sample (required)
        2.  Group: Sample group information (required)
        3+. Remaining N-columns: Optional covariates for correction or 
            additional metadata for eda purposes.
    The contrasts file is a tab-delimited file containing two columns:
        1.  Experimental group of interest (i.e KO)
        2.  Baseline group (i.e WT)
    @Input:
        Groups and contrasts file (gather),
        Raw Counts Matrix
    @Output:
        Differential gene expression results for all comparisons,
        Interactive HTML report with an overview of the results,
        Normalized counts matrix (TMM-scaled, log2CPM counts),
        Normalized batch-corrected counts matrix (if covariates),
        For each comparison: 
         • Complete, unfiltered DGE results (TSV file),
         • Static (PDF, jpeg) and interactive volcano plots (HTML)
         • Linear model coefficient (used to calculate fold-changes)
    """
    input:
        groups=grps_file,
        contrasts=cmps_file,
        counts=join(workpath,degall_dir,"RSEM_genes_expected_counts.tsv"),
    output:
        render_script=join(workpath,"differential_gene_expression",batch_id,"limma","render.R"),
        html_report=join(workpath,"differential_gene_expression",batch_id,"limma","Limma_DGE_Report.html"),
        # Per-comparison deg results
        deg_results=expand(
            join(workpath,"differential_gene_expression",batch_id,"limma","{case}_vs_{control}","limma_voom_deg_{case}_vs_{control}.tsv"),
            zip, case=case_groups, control=ctrl_groups
        ),
    params:
        rname='pl:limma_genes',
        rmd=join(workpath,"workflow","scripts","limma_dge_report.Rmd"),
        outdir=join(workpath,"differential_gene_expression",batch_id,"limma"),
        # Optional covariates to correct for,
        # either a comma delimited list of
        # columns in the groups file or
        # NULL if no covariates
        covariates_option_value = lambda _: "'{0}'".format(
            covariates,
        ) if covariates else 'NULL',
    container: config['images']['dge']
    shell: textwrap.dedent("""
    # Create wrapper script to render
    # the Rmd into a HTML report
    cat << EOF > {output.render_script}
    #!/usr/bin/env Rscript
    library(rmarkdown)
    rmarkdown::render(
        '{params.rmd}',
        output_file = '{output.html_report}',
        params = list(
          raw = '{input.counts}',
          group = '{input.groups}',
          contrasts = '{input.contrasts}', 
          covariates = {params.covariates_option_value},
          prefix = 'Limma_report',
          wdir = '{params.outdir}'
        )
    )
    EOF
    # Unset any R-related env variables
    # that may have been inherited from
    # the user environment
    unset R_LIBS_USER; unset R_LIBS_SITE;
    # Perform differential gene expression analysis
    Rscript {output.render_script}
    """)
