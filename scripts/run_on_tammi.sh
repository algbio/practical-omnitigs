#!/bin/bash

set -e

if [ $# -eq 0 ]; then
    echo "Error: no snakemake target given"
    exit 1
fi

SOURCED_SEBSCHMI_BASHRC="false"
source /home/sebschmi/.profile
cd /home/sebschmi/git/practical-omnitigs

conda activate practical-omnitigs

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

echo "Arguments: $@" >> "$LOGDIR/run_on_tammi.log"
nohup snakemake --profile config/tammi "$@" >> "$LOGDIR/run_on_tammi.log" 2>&1 &

echo "Started snakemake in background with PID $!"

