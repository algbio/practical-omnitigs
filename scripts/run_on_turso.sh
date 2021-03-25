#!/bin/bash

set -e

if [ $# -eq 0 ]; then
    echo "Error: no snakemake target given"
    exit 1
fi

source /home/sebschmi/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs

# Remove erroneous outputs from previous run
echo "Removing erroneous outputs from previous run"
scripts/delete_erroneous_outputs.py

# Create log directory
LOGDIR="logs/$(date +"%FT%X")/"
mkdir -p "$LOGDIR"
LATEST_LOGDIR_SYMLINK="logs/latest"
rm -f "$LATEST_LOGDIR_SYMLINK"
ln -sr "$LOGDIR" "$LATEST_LOGDIR_SYMLINK"

echo "Storing logs in directory $LOGDIR"
echo "Also symlinked as $LATEST_LOGDIR_SYMLINK"
echo "$LOGDIR" > .logdir

echo "Creating jobs"

echo "Arguments: $@" >> "$LOGDIR/run_on_turso.log"
nohup snakemake --profile config/turso "$@" >> "$LOGDIR/run_on_turso.log" 2>&1 &

echo "Started snakemake in background with PID $!"



### Below is just for slurm-based scheduling, which never really worked.

exit 0

nohup snakemake --profile config/turso $@ 2>&1 | tee "$LOGDIR/run_on_turso.log"
rm -f .logdir


echo "Created jobs, logging their properties"

# Log job information
for JOB in $(squeue -o "%.18A" -u sebschmi -M carrington | sed 's/ *//' | grep -E "^[0-9]*$"); do
	scontrol show jobid -M carrington -dd $JOB >> "$LOGDIR/jobs.log"
done
for JOB in $(squeue -o "%.18A" -u sebschmi -M ukko2 | sed 's/ *//' | grep -E "^[0-9]*$"); do
	scontrol show jobid -M ukko2 -dd $JOB >> "$LOGDIR/jobs.log"
done
for JOB in $(squeue -o "%.18A" -u sebschmi -M vorna | sed 's/ *//' | grep -E "^[0-9]*$"); do
	scontrol show jobid -M vorna -dd $JOB >> "$LOGDIR/jobs.log"
done

echo "Successfully created jobs, now releasing them"

squeue -o "%.18A %.18R" -u sebschmi -M carrington | awk '{if ($2 =="(JobHeldUser)"){printf "%s\n", $1}}' | xargs -r scontrol -M carrington release
squeue -o "%.18A %.18R" -u sebschmi -M ukko2 | awk '{if ($2 =="(JobHeldUser)"){printf "%s\n", $1}}' | xargs -r scontrol -M ukko2 release
squeue -o "%.18A %.18R" -u sebschmi -M vorna | awk '{if ($2 =="(JobHeldUser)"){printf "%s\n", $1}}' | xargs -r scontrol -M vorna release

#rm -f .tmpjobids
#squeue -o "%.18A %.18R" -u sebschmi -M carrington | awk '{if ($2 =="(JobHeldUser)"){printf "%s carrington\n", $1}}' >> .tmpjobids
#squeue -o "%.18A %.18R" -u sebschmi -M ukko2 | awk '{if ($2 =="(JobHeldUser)"){printf "%s ukko2\n", $1}}' >> .tmpjobids
#squeue -o "%.18A %.18R" -u sebschmi -M vorna | awk '{if ($2 =="(JobHeldUser)"){printf "%s vorna\n", $1}}' >> .tmpjobids

#sort -nr .tmpjobids | xargs -nr2 sh -c 'scontrol -M $2 release $1' sh
#rm -f .tmpjobids

BROKEN_JOBS=$(squeue -o "%.18A %.18R" -u sebschmi | grep "(DependencyNeverSatisfied)")

if [ -z "$BROKEN_JOBS" ]; then
	echo "Successfully created and released all jobs"
else
	echo "Found broken jobs: $BROKEN_JOBS"
	exit 1
fi
