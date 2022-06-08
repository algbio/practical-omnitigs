#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys

def generate():
    with open(sys.argv[1], 'r') as infile:
        for record in SeqIO.parse(infile, "fasta"):
            last_character = 'N'
            compressed_sequence = []

            for character in record.seq:
                if character != last_character:
                    compressed_sequence.append(character)
                    last_character = character

            record.seq = Seq("".join(compressed_sequence))
            yield record

with open(sys.argv[2], 'w') as outfile:
    SeqIO.write(generate(), outfile, "fasta")