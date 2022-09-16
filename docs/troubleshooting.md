If you are experiencing an issue, please read through this list first before contacting our team.

We have compiled this FAQ from the most common problems. If you are running into an issue that is not on this page, please feel free to [reach out to our team](https://github.com/skchronicles/RNA-seek/issues).

## Job Status

**Q. How do I know if RNA-seek finished running successfully?**

**A.** There are several different ways of checking the status of each job submitted to the cluster.  
Here are a few suggestions:

!!! tldr "Check Job Status"

    === "Biowulf Dashboard"

        You can check the status of Biowulf jobs through the your [user dashboard](https://hpc.nih.gov/dashboard/).

        Each job that RNA-seek submits to the cluster starts with the `pl:` prefix.

    === "Snakemake Log"

        [Snakemake](https://snakemake.readthedocs.io/en/stable/) generates the following file, `Reports/snakemake.log`, in each pipeline's working directory. This file contains information about each job submitted to the job scheduler. If there are no problems, snakemake will report 100% steps done in those last few lines of the file.

        You can take a peek of the end of the file by running the following command:
        ```bash
        tail -n30 Reports/snakemake.log
        ```

        Or more specifically, you can pull out the timestamps of the last few completed jobs like this:
        ```bash
        grep -A 1 done Reports/snakemake.log | tail
        ```

    === "Query Job Scheduler"

        SLURM has built-in commands that allow a user to view the status of jobs submitted to the cluster.

        **Method 1:** To see what jobs you have running, run the following command:
        ```bash
        squeue -u $USER
        ```

        **Method 2** You can also run this alternative command to check the status of your running jobs:
        ```bash
        sjobs
        ```

        Each job that RNA-seek submits to the cluster starts with the `pl:` prefix.


**Q. How do I identify failed jobs?**

**A.** If there are errors, you'll need to identify which jobs failed and check its corresponding SLURM output file.
The SLURM output file may contain a clue as to why the job failed.

!!! tldr "Find Failed Jobs"

    === "SLURM output files"

        Quick and dirty method to search for failed jobs by looking through each job's output file:

        ```bash
        grep -i 'fail' slurmfiles/slurm-*.out
        ```

    === "Snakemake Log"

        [Bash script]( https://github.com/CCBR/Tools/blob/master/Biowulf/get_slurm_file_with_error.sh) identify the SLURM ID of the first failed job and check if the output file exists.


Many failures are caused by filesystem or network issues on Biowulf, and in such cases, simply re-starting the Pipeline should resolve the issue. Snakemake will dynamically determine which steps have been completed, and which steps still need to be run. If you are still running into problems after re-running the pipeline, there may be another issue. If that is the case, please feel free to [contact us](https://github.com/skchronicles/RNA-seek/issues).


**Q. How do I cancel ongoing RNA-seek jobs?**

**A.** Sometimes, you might need to manually stop a RNA-seek run prematurely, perhaps because the run was configured incorrectly or if a job is stalled.  Although the walltime limits will eventually stop the workflow, this can take up to 5 or 10 days depending on the pipeline.

To stop RNA-seek jobs that are currently running, you can follow these options.

!!! tldr "Cancel running jobs"

    === "Master Job"
        You can use the `sjobs` tool [provided by Biowulf](https://hpc.nih.gov/docs/biowulf_tools.html#sjobs) to monitor ongoing jobs.

        Examine the `NAME` column of the `sjobs` output, one of them should match `pl:rna-seek`. This is the "primary" job that orchestrates the submission of child jobs as the pipeline completes. Terminating this job will ensure that the pipeline is cancelled; however, you will likely need to unlock the working directory before re-running rna-seek again. Please see our instructions below in `Error: Directory cannot be locked` for how to unlock a working directory.

        You can [manually cancel](https://hpc.nih.gov/docs/userguide.html#delete) the primary job using `scancel`.

        However, secondary jobs that are already running will continue to completion (or failure).  To stop them immediately, you will need to run `scancel` individually for each secondary job. See the next tab for a bash script that tries to automate this process.

    === "Child Jobs"
        When there are lots of secondary jobs running, or if you have multiple RNA-seek runs ongoing simultaneously, it's not feasible to manually cancel jobs based on the `sjobs` output (see previous tab).

        We provide [a script](https://github.com/CCBR/Tools/blob/master/Biowulf/cancel_snakemake_jobs.sh) that will parse the snakemake log file and cancel all jobs listed within.

        ```bash
        ## Download the script (to the current directory)
        wget https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/cancel_snakemake_jobs.sh

        ## Run the script
        bash cancel_snakemake_jobs.sh /path/to/output/logfiles/snakemake.log
        ```

        The script accepts one argument, which should be the path to the snakemake log file.  This will work for any log output from Snakemake.

        This script will NOT cancel the primary job, which you will still have to identify and cancel manually, as described in the previous tab.

Once you've ensured that all running jobs have been stopped, you need to unlock the working directory (see below), and re-run rna-seek to resume the pipeline.

## Job Errors

**Q. Why am I getting `sbatch: command not found error`?**

**A.** Are you running the `rna-seek` on `helix.nih.gov` by mistake. [Helix](https://hpc.nih.gov/systems/) does not have a job scheduler. One may be able to fire up the singularity module, initial working directory and perform dry-run on `helix`. But to submit jobs, you need to log into `biowulf` using `ssh -Y username@biowulf.nih.gov`.


**Q. Why am I getting a message saying `Error: Directory cannot be locked. ...` when I do the dry-run?**

**A.** This is caused when a run is stopped prematurely, either accidentally or on purpose, or the pipeline is still running in your working directory. Snakemake will lock a working directory to prevent two concurrent pipelines from writing to the same location. This can be remedied easily by running `rna-seek unlock` sub command. Please check to see if the pipeline is still running prior to running the commands below. If you would like to cancel a submitted or running pipeline, please reference the instructions above.

```bash
# Load Dependencies
modue load singularity snakemake

# Unlock the working directory
rna-seek unlock --output /path/to/working/dir
```

**Q. Why am I getting a message saying `MissingInputException in line ...` when I do the dry-run?**

**A.** This error usually occurs when snakemake is terminated ungracefully. Did you forcefully cancel a running pipeline? Or did one of your running pipelines abruptly end? Either way, the solution is straight-forward. Please go to your pipeline's output directory, and rename or delete the following hidden directory: `.snakemake/`. This directory contains metadata pertaining any snakemake runs inside that working directory. Sometimes when a pipeline is pre-maturely or forcefully terminated, a few files related to tracking temp() files are not deleted and snakemake raises a MissingInputException.  

```bash
# Navigate to working directory
cd /path/to/working/dir

# Rename .snakemake directory to something else
# And try re-dry running the pipeline
mv .snakemake .old_snakemake
```

