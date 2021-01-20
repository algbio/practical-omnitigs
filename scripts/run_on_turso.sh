#!/bin/bash

set -e

source /home/sebschmi/.bashrc
cd /proj/sebschmi/git/practical-omnitigs

source activate practical-omnitigs

snakemake --profile config/turso report_all
