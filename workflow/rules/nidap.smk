localrules: nidap
rule nidap:
    input:
        get_nidap_folder_input_files
    output:
        join(workpath,"NIDAP","RNA_Report.html"),
        join(workpath,"NIDAP","multiqc_report.html"),
        join(workpath,"NIDAP","RSEM.genes.FPKM.all_samples.txt"),
        join(workpath,"NIDAP","RSEM.isoforms.FPKM.all_samples.txt"),
    params:
        outdir=join(workpath,"NIDAP")
    shell:"""
set -exo pipefail
if [ -d {params.outdir} ];then rm -rf {params.outdir};fi 
mkdir -p {params.outdir}
cd {params.outdir}
for input in {input};do
    ln $input .
done
"""