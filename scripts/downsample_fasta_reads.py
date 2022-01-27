#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import random

input_file = sys.argv[1]
output_file = sys.argv[2]
factor = float(sys.argv[3])

# count reads
with open(input_file, 'r') as infile:
    count = len(SeqIO.parse(infile, "fasta"))

indices = [i for i in 0:count]
random.Random(3252457).shuffle(indices)
indices = set(indices[0:int(count * factor)])

def generate():
    with open(sys.argv[1], 'r') as infile:
        for index, record in enumerate(SeqIO.parse(infile, "fasta")):
            if index in indices:
                yield record

with open(output_file, 'w') as outfile:
    SeqIO.write(generate(), outfile, "fasta")