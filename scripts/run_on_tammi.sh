#!/bin/bash

set -e

source /home/sebschmi/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs

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
nohup snakemake --profile config/turso $@ >> "$LOGDIR/run_on_turso.log" 2>&1 &

echo "Started snakemake in background with PID $!"

