#!/bin/bash

set -e

source /home/sebschmi/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs

snakemake --profile config/turso $@ | tee run_on_turso.log

squeue -o "%.18A %.18R" -u sebschmi | awk '{if ($2 =="(JobHeldUser)"){print $1}}' | xargs -n 1 scontrol resume

echo "Successfully created and released all jobs"
