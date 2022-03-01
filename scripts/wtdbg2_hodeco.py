#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os
import random

input_file = sys.argv[1]
output_file = sys.argv[2]
temp_file = output_file + ".tmp_ctg"
normal_reads_file = sys.argv[3]
hoco_reads_file = sys.argv[4]

# load reads
print(f"Loading normal reads from: {normal_reads_file}", flush = True)
normal_reads = {}
with open(normal_reads_file, 'r') as infile:
    for record in SeqIO.parse(infile, "fasta"):
        normal_reads[record.id] = record.seq

print(f"Loading hoco reads from: {hoco_reads_file}", flush = True)
hoco_reads = {}
with open(hoco_reads_file, 'r') as infile:
    for record in SeqIO.parse(infile, "fasta"):
        hoco_reads[record.id] = record.seq

def transform_indices_to_hoco_forwards(read_id, read_offset, read_limit, modified_offset_shifts, lay_counts):
    normal_sequence = normal_reads[read_id]
    hoco_sequence = hoco_reads[read_id]
    shifted_read_offset = 0
    shifted_read_limit = 0

    for i in range(0, read_offset):
        assert hoco_sequence[i] == normal_sequence[shifted_read_offset]
        assert shifted_read_offset < len(normal_sequence)
        while shifted_read_offset < len(normal_sequence) and hoco_sequence[i] == normal_sequence[shifted_read_offset]:
            shifted_read_offset += 1

    shifted_read_limit = shifted_read_offset
    for i in range(read_offset, read_limit):
        hoco_align_index = i - read_offset
        while len(modified_offset_shifts) <= hoco_align_index:
            modified_offset_shifts.append(0)
            lay_counts.append(0)
        modified_offset_shifts[hoco_align_index] += (shifted_read_limit - shifted_read_offset) - hoco_align_index
        lay_counts[hoco_align_index] += 1

        assert hoco_sequence[i] == normal_sequence[shifted_read_limit]
        assert shifted_read_limit < len(normal_sequence)
        while shifted_read_limit < len(normal_sequence) and hoco_sequence[i] == normal_sequence[shifted_read_limit]:
            shifted_read_limit += 1

    return shifted_read_offset, shifted_read_limit

def transform_indices_to_hoco_backwards(read_id, read_offset, read_limit, modified_offset_shifts, lay_counts):
    normal_sequence = normal_reads[read_id]
    hoco_sequence = hoco_reads[read_id]
    shifted_read_offset = len(normal_sequence) - 1
    shifted_read_limit = len(normal_sequence) - 1

    for i in reversed(range(read_limit, len(hoco_sequence))):
        assert hoco_sequence[i] == normal_sequence[shifted_read_limit]
        assert shifted_read_limit >= 0
        while hoco_sequence[i] == normal_sequence[shifted_read_limit] and shifted_read_limit > 0:
            shifted_read_limit -= 1

    shifted_read_offset = shifted_read_limit
    shifted_read_limit += 1
    for i in reversed(range(read_offset, read_limit)):
        hoco_align_index = read_limit - i - 1
        while len(modified_offset_shifts) <= hoco_align_index:
            modified_offset_shifts.append(0)
            lay_counts.append(0)
        modified_offset_shifts[hoco_align_index] += (shifted_read_limit - shifted_read_offset) - hoco_align_index
        lay_counts[hoco_align_index] += 1

        assert hoco_sequence[i] == normal_sequence[shifted_read_offset]
        assert shifted_read_offset >= 0
        while hoco_sequence[i] == normal_sequence[shifted_read_offset] and shifted_read_offset > 0:
            shifted_read_offset -= 1
    shifted_read_offset += 1

    return shifted_read_offset, shifted_read_limit

modified_offset_shifts = []
lay_counts = []

# process contigs
print("Decompressing contigs", flush = True)
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        line = line.strip()
        if line.startswith('>'):
            if len(modified_offset_shifts) > 0:
                assert len(modified_offset_shifts) == len(lay_counts)
                current_offset_shift += int(modified_offset_shifts[len(modified_offset_shifts) - 1] / lay_counts[len(lay_counts) - 1])
                modified_offset_shifts = []
                lay_counts = []

                outfile.write(f">{contig_id}\tnodes={contig_nodes}\tlen={contig_len + current_offset_shift}\n")
                tmpfile.close()
                tmpfile = open(temp_file, 'r')
                while True:
                    data = tmpfile.read(1024 * 1024)
                    if data:
                        while len(data) > 0:
                            written = outfile.write(data)
                            data = data[written:]
                    else:
                        break
                tmpfile.close()
                os.remove(temp_file)

            tmpfile = open(temp_file, 'w')
            line = line[1:]
            line = list(line.split(' '))
            contig_id = line[0]
            contig_nodes = int(line[1][6:])
            contig_len = int(line[2][4:])

            print(f"Processing contig {contig_id} of hoco len: {contig_len}", flush = True)
            current_offset_shift = 0
            modified_offset_shifts = []
            lay_counts = []
            old_input_offset = None
            input_offset = None
        elif line.startswith('E'):
            line = list(line.split('\t'))
            old_input_offset = input_offset
            input_offset = int(line[1])
            n1 = line[2]
            n1_direction = line[3]
            n2 = line[4]
            n2_direction = line[5]

            if len(modified_offset_shifts) > 0:
                assert len(modified_offset_shifts) == len(lay_counts)
                index = min(len(modified_offset_shifts) - 1, input_offset - old_input_offset)
                current_offset_shift += int(modified_offset_shifts[index] / lay_counts[index])
                modified_offset_shifts = []
                lay_counts = []

            tmpfile.write(f"E\t{input_offset + current_offset_shift}\t{n1}\t{n1_direction}\t{n2}\t{n2_direction}\n")
        elif line.startswith('S'):
            line = list(line.split('\t'))
            read_id = line[1]
            read_direction = line[2]
            read_offset = int(line[3])
            read_len = int(line[4])
            read_limit = read_offset + read_len
            read_sequence = line[5]

            # transform hoco indices to normal indices
            if read_direction == '+':
                shifted_read_offset, shifted_read_limit = transform_indices_to_hoco_forwards(read_id, read_offset, read_limit, modified_offset_shifts, lay_counts)
            elif read_direction == '-':
                shifted_read_offset, shifted_read_limit = transform_indices_to_hoco_backwards(read_id, read_offset, read_limit, modified_offset_shifts, lay_counts)
            else:
                raise Exception(f"Unknown read direction: {read_direction}")

            shifted_read_len = shifted_read_limit - shifted_read_offset
            shifted_read_sequence = normal_reads[read_id][shifted_read_offset:shifted_read_limit]
            if read_direction == '-':
                shifted_read_sequence = shifted_read_sequence.reverse_complement()

            tmpfile.write(f"S\t{read_id}\t{read_direction}\t{shifted_read_offset + current_offset_shift}\t{shifted_read_len}\t{shifted_read_sequence}\n")
        else:
            raise Exception(f"Unknown line start: {line}")

    if len(modified_offset_shifts) > 0:
        assert len(modified_offset_shifts) == len(lay_counts)
        current_offset_shift += int(modified_offset_shifts[len(modified_offset_shifts) - 1] / lay_counts[len(lay_counts) - 1])
        modified_offset_shifts = []
        lay_counts = []

        outfile.write(f">{contig_id}\tnodes={contig_nodes}\tlen={contig_len + current_offset_shift}\n")
        tmpfile.close()
        tmpfile = open(temp_file, 'r')
        while True:
            data = tmpfile.read(1024 * 1024)
            if data:
                while len(data) > 0:
                    written = outfile.write(data)
                    data = data[written:]
            else:
                break
        tmpfile.close()
        os.remove(temp_file)

print("Done", flush = True)