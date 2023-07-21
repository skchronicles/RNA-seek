#!/usr/bin/env bash
set -euo pipefail

# USAGE: ./do_create_rMATS_sample_sheet.sh --skip-index

# Functions
function err() { cat <<< "$@" 1>&2; }
function usage() { cat << EOF
do_run_rMATS: Run rMATS with an RNA-seek output directory.
USAGE:
  $ ./do_run_rMATS [-h] \\
        [--skip-index]
SYNOPSIS:
  Convience script to run rMATS with an RNA-seek output directory. A 
  user just needs to create a 'groups.tab' and 'contrasts.tab' file in the
  RNA-seek output directory of interest. This script will run rMATS Turbo
  for each comparsion defined in 'contrasts.tab' file. The '--skip-index'
  option can be provided if this script has already been run to generate 
  an STAR index for rMATS.
OPTIONS:
  -s, --skip-index [Type: Bool] Skip building STAR's Index. 
                                 WARNING: there be dragens here! This option 
                                 should only be provided if ./do_run_rMATS has
                                 already built an index in the specified output
                                 directory in the past. You may want to provide 
                                 this option if rMATS failed to run for some
                                 other weird issue AND you do not want to waste
                                 time rebuilding the index. Do not provide this
                                 option if you do not know what you are doing! 
                                 Running rMATS without its index for STAR will 
                                 cause rMATS to fail. 
  -h, --help     [Type: Bool]  Displays usage and help information.
Example:
  $ cd /path/to/RNA-seek/output/directory
  $ nano groups.tab    # create a groups.tab file similar to Pipeliner 
  $ nano contrasts.tab # create a contrasts.tab file similar to Pipeliner
  $ ./do_run_rMATS
Version:
  0.2.0
EOF
}


function create_groups(){
    # Creates required sample sheet to detect differential 
    # AS events with rMATS turbo. 
    
    # Create sample sheet from groups.tab and contrasts.tab
    while read g1 g2; do 
        # Get list of samples for the first group
        s1=$(awk \
                -F '\t' \
                -v group="$g1" \
                -v wd="$PWD" \
                '$2==group {print wd"/"$1".R1.fastq.gz:"wd"/"$1".R2.fastq.gz"}' \
            groups.tab | tr '\n' ',');
        # Get list of samples for the second group
        s2=$(awk \
                -F '\t' \
                -v group="$g2" \
                -v wd="$PWD" \
                '$2==group {print wd"/"$1".R1.fastq.gz:"wd"/"$1".R2.fastq.gz"}' \
            groups.tab | tr '\n' ','); 
        # Create sample sheet for the first group
        echo -e  "${s1%,}" > "rMATS/${g1}_v_${g2}/s1.txt"; 
        # Create sample sheet for the second group
        echo -e  "${s2%,}" > "rMATS/${g1}_v_${g2}/s2.txt"; 
    done < contrasts.tab
}


function initalize(){
    # Creates output directory heirarchy
    # $1 = output directory
    local wd="$1"

    ( # Initialize a base output directory structure 
        cd ${wd};
        mkdir -p "${wd}/rMATS"
        while read g1 g2; do 
            # Create an output directory for each contrast
            mkdir -p "${wd}/rMATS/${g1}_v_${g2}/";
        done < contrasts.tab
    )
}


function _get_read_length(){
    # Get max read length for rMATS model
    # $1 = MultiQC matrix with read lengths 
    local metadata="$1"
    cut -f6 "$metadata" \
        | tail -n+2 \
        | awk -F '-' '{print $NF}' \
        | sort -k1,1nr \
        | head -1
}


function build_star_index(){
    # Submits a job to builds an index for STAR
    # @RETURNS sbatch JobID of build Index job

    local read_length
    local gtf
    local genome
    local star_version="2.7.6a"

    # Required for parsing config file
    module load jq > /dev/null 2>&1

    # Get maximum read length 
    # and create output directory for new index
    read_length=$(_get_read_length "Reports/multiqc_matrix.tsv")
    mkdir -p rMATS/STAR/2.7.6a/genes-${read_length}

    # Get references to build STAR index
    gtf=$(jq .references.rnaseq.GTFFILE config.json)
    genome=$(jq .references.rnaseq.GENOME config.json)
    rl=$((${read_length}-1))
    
    # Create sbacth script to build index
    cat << EOF > build_star_index_submit.sh
#!/usr/bin/env bash
#SBATCH --cpus-per-task=32 
#SBATCH --mem=64g
#SBATCH --gres=lscratch:250 
#SBATCH --time=8:00:00 
#SBATCH --parsable 
#SBATCH -J "STAR_INDEX" 
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
module load STAR/${star_version}

# Builds STAR Index to align reads against reference genome for a defined
# read length. Optimal readlength for sjdbOverhang = max(ReadLength) - 1.
# For most applications, a readlength of 100 works well.
STAR \\
	--runThreadN 32 \\
	--runMode genomeGenerate \\
	--genomeDir rMATS/STAR/2.7.6a/genes-${read_length} \\
	--genomeFastaFiles ${genome} \\
	--sjdbGTFfile ${gtf} \\
	--sjdbOverhang ${rl} \\
	--outFileNamePrefix rMATS/STAR/2.7.6a/build_${read_length}_ \\
	--outTmpDir /lscratch/\${SLURM_JOB_ID}/tmp_${read_length}
EOF
    err "Submitting job to build STAR Index"
    chmod +x build_star_index_submit.sh
    sbatch build_star_index_submit.sh
}


function do_run_rMATS(){
    # Submits job to run rMATS Turbo a contrast
    # $1 = Group 1
    # $2 = Group 2
    # $3 = STAR Index 
    # $4 = Output directory
    # $5 = Job dependency (slurm job id of build index)
    local gtf
    local read_length
    local g1="${1}" 
    local g2="${2}" 
    local star_index="${3}"
    local outdir="${4}"
    local dependency="${5:-}"
    local star_version="2.7.6a"

    # Required for parsing config file
    module load jq > /dev/null 2>&1
    gtf=$(jq .references.rnaseq.GTFFILE config.json)
    # Get maximum read length 
    read_length=$(_get_read_length "Reports/multiqc_matrix.tsv")

    # Create sbacth script to run rMATS
    cat << EOF > run_rmats_${g1}_${g2}.sh
#!/usr/bin/env bash
#SBATCH --cpus-per-task=32 
#SBATCH --mem=64g
#SBATCH --gres=lscratch:250 
#SBATCH --time=8:00:00 
#SBATCH --parsable 
#SBATCH -J "rMATS_${g1}_${g2}" 
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
module load rMATS/4.1.2 
module load STAR/${star_version}

# Run rMATS for a given contrast 
# starting from untrimmed FastQ files
rmats.py \\
    --s1 "${PWD}/rMATS/${g1}_v_${g2}/s1.txt" \\
    --s2 "${PWD}/rMATS/${g1}_v_${g2}/s2.txt" \\
    --od "${outdir}" \\
    --gtf ${gtf} \\
    --readLength ${read_length} \\
    --bi "${star_index}" \\
    --nthread 32 \\
    -t "paired" \\
    --tmp /lscratch/\${SLURM_JOB_ID}/
EOF
    
    chmod +x run_rmats_${g1}_${g2}.sh
    if [ -z "${dependency}" ]; then
        # No job dependency
        sbatch run_rmats_${g1}_${g2}.sh
    else
        # Submit with STAR index build dependency
        sbatch --dependency=afterok:${dependency} run_rmats_${g1}_${g2}.sh
    fi
}

function main(){
    # Pseudo main method
    # Change directory to script's working directory 
    # (i.e pipeline output directory)
    cd "$(dirname "${BASH_SOURCE[0]}")"

    # Parser command line arguments 
    # Associative array to store parsed args
    declare -Ag args
    args["skip_index"]=false # default behavior does build index!
    while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
          -h  | --help) usage && exit 0;;
          -b | --skip-index) args["skip_index"]=true; shift;;
          -*  | --*) err "Error: Failed to parse unsupported argument: '${key}'."; usage && exit 1;;
          *) err "Error: Failed to parse unrecognized argument: '${key}'. Do any of your inputs have spaces?"; usage && exit 1;;
        esac
    done

    # Step 1. Create a base output directory hierarchy
    initalize "${PWD}"
    
    # Step 2. Create samples sheet for each contrast
    # from a groups.tab file and contrasts.tab file
    create_groups 
    
    # Step 3. Build STAR Index 
    if [ "${args["skip_index"]}" = false ]; then
        build_job_id=$(build_star_index)
    fi

    # Get index for the correct read length 
    read_length_index=$(_get_read_length "Reports/multiqc_matrix.tsv" \
                        | awk -v wd="$PWD" \
                            '{print wd"/rMATS/STAR/2.7.6a/genes-"$1"/"}')
    while read g1 g2; do 
        # Run rMATS for each constrast
        # $1 = Group 1
        # $2 = Group 2
        # $3 = STAR Index 
        # $4 = Output directory
        # $5 = Job dependency (slurm job id of build index)
        do_run_rMATS "${g1}" "${g2}" "${read_length_index}" "${PWD}/rMATS/${g1}_v_${g2}/" "${build_job_id:-}"

    done < contrasts.tab

}


main "$@"
