#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys

def generate():
    with open(sys.argv[1], 'r') as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if not record.id.startswith("CM"):
                continue
            yield record

with open(sys.argv[2], 'w') as outfile:
    SeqIO.write(generate(), outfile, "fasta")
