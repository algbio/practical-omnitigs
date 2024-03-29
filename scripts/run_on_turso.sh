#!/bin/bash

set +e

echo "Checking for running snakemake"
ps ux | grep -F "practical-omnitigs/bin/snakemake " | grep -v grep
if [ "${PIPESTATUS[2]}" -ne "0" ]; then
    echo "No snakemake found"
else
    echo "Snakemake is already running!"
    exit 1
fi

echo "Checking if there are any existing slurm jobs"
if [ -z "$(squeue -o '%A %.28R %j' -u $(whoami) -M all | grep 'practical-omnitigs')" ]; then
    echo "No slurm jobs found"
else
    echo "Found slurm jobs"
    squeue -o '%A %.28R %j' -u $(whoami) -M all | grep 'practical-omnitigs'
    exit 1
fi

set -e

if [ $# -eq 0 ]; then
    echo "Error: no snakemake target given"
    exit 1
fi

source $HOME/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs

# Remove erroneous outputs from previous run
echo "Removing erroneous outputs from previous run"
scripts/delete_erroneous_outputs.py

# Create log directory
LOGDIR="logs/$(date +"%FT%X")/"
mkdir -p "$LOGDIR"
LOGDIR=$(realpath $LOGDIR)
LATEST_LOGDIR_SYMLINK="logs/latest"
rm -f "$LATEST_LOGDIR_SYMLINK"
ln -sr "$LOGDIR" "$LATEST_LOGDIR_SYMLINK"

echo "Storing logs in directory $LOGDIR"
echo "Also symlinked as $LATEST_LOGDIR_SYMLINK"
echo "$LOGDIR" > .logdir

# Clean shadow directory
# echo "Cleaning shadow directory"
# snakemake --profile config/turso --cleanup-shadow >> "$LOGDIR/run_on_turso.log" 2>&1

export TMPDIR="/wrk-vakka/users/sebschmi/practical-omnitigs/tmp"
export TEMP=$TMPDIR
export TMP=$TMPDIR
mkdir -p "$TMPDIR"

# Build environments (due to a bug in snakemake, this needs to be done in a separate step for now)
# echo "Building conda environments"
# snakemake --use-conda --cores 1 --conda-create-envs-only "$@" 2>&1 | tee "$LOGDIR/run_on_turso.log"

# Execute pipeline
echo "Creating jobs"

echo "Arguments: $@" >> "$LOGDIR/run_on_turso.log"
nohup snakemake --profile config/turso "$@" >> "$LOGDIR/run_on_turso.log" 2>&1 &

echo "Started snakemake in background with PID $!"
