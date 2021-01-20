#!/bin/bash

set -e

source /home/sebschmi/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs

snakemake --profile config/turso $@ | tee run_on_turso.log

squeue -o "%.18A %.18R" -u sebschmi -M carrington | awk '{if ($2 =="(JobHeldUser)"){print $1}}' | xargs -n 1 -r scontrol -M carrington release
squeue -o "%.18A %.18R" -u sebschmi -M ukko2 | awk '{if ($2 =="(JobHeldUser)"){print $1}}' | xargs -n 1 -r scontrol -M ukko2 release
squeue -o "%.18A %.18R" -u sebschmi -M vorna | awk '{if ($2 =="(JobHeldUser)"){print $1}}' | xargs -n 1 -r scontrol -M vorna release

echo "Successfully created and released all jobs"
