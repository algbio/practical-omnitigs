#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys

min_read_len = None
with open(sys.argv[1], 'r') as infile:
    for record in SeqIO.parse(infile, "fasta"):
        read_len = len(record.seq)
        if min_read_len is None:
            min_read_len = read_len
        else:
            min_read_len = min(min_read_len, read_len)

print(f"Minimum read length: {min_read_len}")
