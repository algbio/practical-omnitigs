#!/usr/bin/python3

import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    sys.exit("Need exactly 2 parameters: <input-file> <output-file>")

input_file = sys.argv[1]
output_file = sys.argv[2]

def read_fasta(input_file):
    id = 0
    for record in SeqIO.parse(input_file, "fasta"):
        record.id = str(id) + "_" + record.id
        id += 1
        yield record

SeqIO.write(read_fasta(input_file), output_file, "fasta")