#!/bin/bash
# helper script that parses slurm output for the job ID,
# and feeds it to back to snakemake/slurm for dependencies.
# This is required when you want to use the snakemake --immediate-submit option

if [[ "Submitted batch job" =~ "$@" ]]; then
  echo -n ""
else
  deplist=$(grep -Eo '[0-9]{5,10}' <<< "$@" | tr '\n' ':' | sed 's/.$//') # | sed 's/,/,afterok:/g')
  #echo "Creating job with deplist $deplist" >> logs/deplists.log
  echo -n "--dependency=afterok:$deplist"
fi
