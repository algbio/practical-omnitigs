#!/bin/bash

set -e

source /home/sebschmi/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs
mkdir -p logs

echo "Creating jobs"

snakemake --profile config/turso $@ | tee "logs/run_on_turso.log"

echo "Created jobs, logging their properties"

# Log job information
rm -f logs/jobs.log
for JOB in $(squeue -o "%.18A" -u sebschmi -M carrington | sed 's/ *//' | grep -E "^[0-9]*$"); do
	scontrol show jobid -M carrington -dd $JOB >> logs/jobs.log
done
for JOB in $(squeue -o "%.18A" -u sebschmi -M ukko2 | sed 's/ *//' | grep -E "^[0-9]*$"); do
	scontrol show jobid -M ukko2 -dd $JOB >> logs/jobs.log
done
for JOB in $(squeue -o "%.18A" -u sebschmi -M vorna | sed 's/ *//' | grep -E "^[0-9]*$"); do
	scontrol show jobid -M vorna -dd $JOB >> logs/jobs.log
done

echo "Successfully created jobs, now releasing them"

rm -f .tmpjobids
squeue -o "%.18A %.18R" -u sebschmi -M carrington | awk '{if ($2 =="(JobHeldUser)"){printf "%s carrington\n", $1}}' >> .tmpjobids
squeue -o "%.18A %.18R" -u sebschmi -M ukko2 | awk '{if ($2 =="(JobHeldUser)"){printf "%s ukko2\n", $1}}' >> .tmpjobids
squeue -o "%.18A %.18R" -u sebschmi -M vorna | awk '{if ($2 =="(JobHeldUser)"){printf "%s vorna\n", $1}}' >> .tmpjobids

sort -nr .tmpjobids | xargs -n2 sh -c 'scontrol -M $2 release $1' sh
rm -f .tmpjobids

BROKEN_JOBS=$(squeue -o "%.18A %.18R" -u sebschmi | grep "(DependencyNeverSatisfied)")

if [ -z "$BROKEN_JOBS" ]; then
	echo "Successfully created and released all jobs"
else
	echo "Found broken jobs: $BROKEN_JOBS"
	exit 1
fi
