#!/usr/bin/env python3

import argparse
import sys
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Simulate perfect reads from a reference")
parser.add_argument("--reference", type=str, required=True, help="Path to the reference file")
parser.add_argument("--output", type=str, required=True, help="Path to the output reads file")
parser.add_argument("--read-length-interval", type=int, nargs=2, default=[15000, 25000], help="Min and max read length")
parser.add_argument("--distribution", type=str, default="cut", help="Sampling method for the reads")
parser.add_argument("--coverage", type=float, default=20.0, help="Target read coverage")

args = parser.parse_args()

reference_path = args.reference
output_path = args.output
read_length_interval = args.read_length_interval
distribution = args.distribution
coverage = args.coverage

def clean_reference(reference):
	for i in range(len(reference)):
		if reference[i] == ord('R'):
			reference[i] = "AG"[random.randint(0, 1)]
		elif reference[i] == ord('Y'):
			reference[i] = "CT"[random.randint(0, 1)]
		elif reference[i] == ord('K'):
			reference[i] = "GT"[random.randint(0, 1)]
		elif reference[i] == ord('M'):
			reference[i] = "AC"[random.randint(0, 1)]
		elif reference[i] == ord('S'):
			reference[i] = "CG"[random.randint(0, 1)]
		elif reference[i] == ord('W'):
			reference[i] = "AT"[random.randint(0, 1)]
		elif reference[i] == ord('B'):
			reference[i] = "CGT"[random.randint(0, 2)]
		elif reference[i] == ord('D'):
			reference[i] = "AGT"[random.randint(0, 2)]
		elif reference[i] == ord('H'):
			reference[i] = "ACT"[random.randint(0, 2)]
		elif reference[i] == ord('V'):
			reference[i] = "ACG"[random.randint(0, 2)]
	return reference

def break_reference(reference):
	offset = 0
	limit = 0

	while limit < len(reference):
		if reference[limit] in b"ACGT":
			pass
		elif reference[limit] == ord('a'):
			reference[limit] = ord('A')
		elif reference[limit] == ord('c'):
			reference[limit] = ord('C')
		elif reference[limit] == ord('g'):
			reference[limit] = ord('G')
		elif reference[limit] == ord('t'):
			reference[limit] = ord('T')
		elif reference[limit] == ord('N') or reference[limit] == ord('-'):
			if limit - offset >= read_length_interval[0]:
				yield offset, reference[offset:limit]
			offset = limit + 1
		else:
			sys.exit("Found unknown character in reference: {}".format(str(reference[limit:limit + 1], "ASCII")))

		limit += 1

	if limit - offset >= read_length_interval[0]:
		yield offset, reference[offset:limit]

def reverse_complement_sequence(sequence):
	reverse = bytearray(sequence)
	for i, c in enumerate(sequence):
		if c == ord('A'):
			rc = ord('T')
		elif c == ord('T'):
			rc = ord('A')
		elif c == ord('C'):
			rc = ord('G')
		elif c == ord('G'):
			rc = ord('C')
		else:
			sys.exit("Found unknown character in sequence: {}".format(str(b"" + c, "ASCII")))
		reverse[len(sequence) - i - 1] = rc
	return reverse

def simulate_cut_reads():
	read_id = 0
	for reference in SeqIO.parse(reference_path, "fasta"):
		reference = clean_reference(bytearray(str(reference.seq).encode("ASCII")))

		for seqence_offset, sequence in break_reference(reference):
			print("Got subsequence of length {}".format(len(sequence)))

			for repetition in range(int(coverage)):
				offset = random.randint(-read_length_interval[0] + 1, 0)

				while offset < len(sequence):
					length = random.randint(read_length_interval[0], read_length_interval[1])
					limit = offset + length
					offset = max(0, offset)
					limit = min(len(sequence), limit)
					if limit - offset >= read_length_interval[0] and limit - offset <= read_length_interval[1]:
						read_seq = sequence[offset:limit]
						if random.randint(0, 1) == 0:
							read_name = str(read_id) + "_" + str(reference.id) + "_" + str(seqence_offset + offset) + "_" + str(seqence_offset + limit)
						else:
							read_seq = reverse_complement_sequence(read_seq)
							read_name = str(read_id) + "_" + str(reference.id) + "_" + str(seqence_offset + limit) + "_" + str(seqence_offset + offset)

						yield SeqRecord(Seq(str(read_seq, "ASCII")), read_name, "", "simulated with cut distribution")
					offset = limit
					read_id += 1

def simulate_reads():
	if distribution == "cut":
		return simulate_cut_reads()
	else:
		sys.exit("Unknown distribution: {}".format(distribution))

SeqIO.write(simulate_reads(), output_path, "fasta")