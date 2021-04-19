#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import logging
import pafpy
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description = "A tool to break contigs at likely breakpoints using long reads.")
parser.add_argument("--input-contigs", metavar = "FASTA_FILE", type = str, required = True, help = "The contigs to break.")
parser.add_argument("--input-reads", metavar = "FASTA_FILE", type = str, required = True, help = "The reads to use for contig breaking.")
parser.add_argument("--output-contigs", metavar = "FASTA_FILE", type = str, required = True, help = "The broken contigs.")
parser.add_argument("--lazy-minimap2", action = "store_true", help = "Do not run minimap2 if its output files already exist.")
parser.add_argument("--blocking-size", metavar = "BLOCKING_SIZE", type = int, default = 32, help = "Process the contigs in blocks of <BLOCKING_SIZE> bases.")
parser.add_argument("--min-surrounding-aligned-bases", metavar = "BASE_COUNT", type = int, default = 500, help = "A read is considered evidence for the correctness of a block if it is surrounded by <BASE_COUNT> aligned bases to the left and right.")
parser.add_argument("--contig-ends-grace-len", metavar = "BASE_COUNT", type = int, default = 700, help = "The first and last <BASE_COUNT> bases of a contig are never considered as breakpoints. Note that it does not make sense to set this value lower than --min-surrounding-aligned-bases.")
parser.add_argument("--min-block-evidence", metavar = "ALIGN_COUNT", type = int, default = 2, help = "A block is considered not a breakpoint if it has at least <ALIGN_COUNT> alignments of evidence.")
parser.add_argument("--min-broken-contig-len", metavar = "BASE_COUNT", type = int, default = 1000, help = "A broken contig is only output if it has at least <BASE_COUNT> bases.")
parser.add_argument("--proper-alignment-contig-ends", metavar = "BASE_COUNT", type = int, default = 10, help = "A read that aligns at most <BASE_COUNT> away from the start or end of a contig is not improperly aligned at that start or end.")
parser.add_argument("--proper-alignment-read-ends", metavar = "BASE_COUNT", type = int, default = 10, help = "A read that aligns at most <BASE_COUNT> away from its start or end is not improperly aligned at that start or end.")
parser.add_argument("--max-ignore-terminating-breakpoints-size", metavar = "BASE_COUNT", type = int, default = 1000, help = "A breakpoint that is at the start or end of a contig and is at most <BASE_COUNT> bases long is ignored.")
parser.add_argument("--threads", metavar = "THREADS", type = int, default = 3, help = "Use <THREADS> threads for minimap2.")
parser.add_argument("--compare-breakpoints", metavar = "JSON_FILE", type = str, help = "Compare with the breakpoints given in <JSON_FILE>.")
parser.add_argument("--algorithm", metavar = "ALGORITHM", type = str, default = "any-inner-alignment", choices = ["any-inner-alignment", "only-proper-alignment"], help = "The algorithm to use for breakpoint detection.")

args = parser.parse_args()

### Set up logging ###

logging_format = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format = logging_format)
logger = logging.getLogger("contigbreaker")
logger.setLevel(logging.INFO)

### Verify arguments ###

if args.min_surrounding_aligned_bases > args.contig_ends_grace_len:
	logger.error("--min-surrounding-aligned-bases is larger than --contig-ends-grace-len. This will result in meaningless breakpoints at the start and end of contigs. Aborting.")
	sys.exit(1)

### Create working directory ###

working_directory = args.output_contigs + ".contigbreakertmp/"
logger.info("Setting up working directory '%s'", working_directory)
if not os.path.isdir(working_directory):
	os.mkdir(working_directory)

### Run minimap2 ###

minimap2_paf_file = working_directory + "alignments.paf"

if os.path.isfile(minimap2_paf_file) and args.lazy_minimap2:
	logger.info("Not running minimap2 because --lazy-minimap2 was specified and the output file already exists: '%s'", minimap2_paf_file)
else:
	if args.lazy_minimap2:
		logger.info("Running minimap2 even though --lazy-minimap2 was specified, because the output file is missing: '%s'", minimap2_paf_file)
	else:
		logger.info("Running minimap2")

	try:
		subprocess.run(["minimap2", "-x", "map-pb", "-O", "1,3", "-p", "0.1", "-N", "10000", "-t", str(args.threads), "-o", minimap2_paf_file, args.input_contigs, args.input_reads], check = True)
	except Exception as e:
		logger.error("Error running minimap2")
		print(e)
		sys.exit(1)

### Load contigs ###

logger.info("Loading contigs from '%s'", args.input_contigs)
contigs = {}

for record in SeqIO.parse(args.input_contigs, "fasta"):
	contigs[record.id] = record

### Load read lenghts ###

logger.info("Loading read lengths from '%s'", args.input_reads)
read_lengths = {}

for record in SeqIO.parse(args.input_reads, "fasta"):
	read_lengths[record.id] = len(record.seq)

### Read paf file ###

logger.info("Reading minimap2 paf file '%s'", minimap2_paf_file)
contig_alignments = {}

with pafpy.PafFile(minimap2_paf_file) as paf_file:
	for record in paf_file:
		contig_alignments.setdefault(record.tname, []).append(record)

### Compute block qualities ###

def div_int_ceil(dividend, divisor):
	return (dividend - 1) // divisor + 1

logger.info("Computing block qualities")

unmapped_reads = 0
unmapped_contigs = set()
total_evidence = 0
total_blocks = 0

contig_block_evidence = {}
for contig_name in contigs.keys():
	contig_len = len(contigs[contig_name].seq)
	block_evidence = contig_block_evidence.setdefault(contig_name, [0.0] * (contig_len // args.blocking_size))

	if contig_name not in contig_alignments:
		unmapped_contigs.add(contig_name)
		continue

	for contig_alignment in contig_alignments[contig_name]:
		if contig_alignment.is_unmapped():
			unmapped_reads += 1
			continue

		if contig_alignment.is_inversion():
			logger.warn("Found inversion, do not know how to handle this: read '%s'", contig_alignment.qname)
			continue

		if contig_alignment.tstart > contig_alignment.tend:
			logger.warn("Found reverse alignment on target, do not know how to handle this: read '%s'", contig_alignment.qname)
			continue

		if args.algorithm == "any-inner-alignment":
			first_block = div_int_ceil(contig_alignment.tstart + args.min_surrounding_aligned_bases, args.blocking_size)
			# The block after the last block
			limit_block = (contig_alignment.tend + 1 - args.min_surrounding_aligned_bases) // args.blocking_size

			for i in range(first_block, limit_block):
				block_evidence[i] += 1
				total_evidence += 1
		elif args.algorithm == "only-proper-alignment":
			read_len = read_lengths[contig_alignment.qname]
			if contig_alignment.tstart > args.proper_alignment_contig_ends and contig_alignment.qstart > args.proper_alignment_read_ends:
				# Alignment does not start at the start of the read
				continue

			if contig_alignment.tend < contig_len - args.proper_alignment_contig_ends and contig_alignment.qend < read_len - args.proper_alignment_read_ends:
				# Alignment does not end at the end of the read
				continue

			first_block = div_int_ceil(contig_alignment.tstart + args.min_surrounding_aligned_bases, args.blocking_size)
			# The block after the last block
			limit_block = (contig_alignment.tend + 1 - args.min_surrounding_aligned_bases) // args.blocking_size

			for i in range(first_block, limit_block):
				block_evidence[i] += 1
				total_evidence += 1

### Compute breakpoints ###

logger.info("Computing breakpoints")
contig_breakpoints = {}
total_breakpoints = 0
non_terminating_breakpoints = 0

for contig_name, block_evidence in contig_block_evidence.items():
	breakpoints = contig_breakpoints.setdefault(contig_name, [])
	contig_len = len(contigs[contig_name].seq)

	first_block = div_int_ceil(args.contig_ends_grace_len, args.blocking_size)
	last_block = (contig_len - args.contig_ends_grace_len) // args.blocking_size
	last_existing_block = contig_len // args.blocking_size
	total_blocks += last_existing_block + 1
	last_breakpoint_start = None

	for i in range(first_block, last_block):
		if block_evidence[i] < args.min_block_evidence:
			if last_breakpoint_start is None:
				last_breakpoint_start = i if i != first_block else 0
		else:
			if last_breakpoint_start is not None:
				if last_breakpoint_start != 0 or i * args.blocking_size > args.max_ignore_terminating_breakpoints_size:
					breakpoints.append((last_breakpoint_start, i - 1))
					total_breakpoints += 1
					if last_breakpoint_start != 0:
						non_terminating_breakpoints += 1
				last_breakpoint_start = None

	if last_breakpoint_start is not None:
		if (last_existing_block - last_breakpoint_start) * args.blocking_size > args.max_ignore_terminating_breakpoints_size:
			breakpoints.append((last_breakpoint_start, last_existing_block))
			total_breakpoints += 1
		last_breakpoint_start = None

	if len(breakpoints) > 0:
		logger.info("Contig '%s' breakpoints: %s", contig_name, ", ".join([str(breakpoint) for breakpoint in breakpoints]))

### Compare breakpoints ###

if args.compare_breakpoints is not None:
	logger.info("Comparing breakpoints to real breakpoints given as '%s'", args.compare_breakpoints)	
	total_correct_breakpoints = 0
	total_wrong_breakpoints = 0
	total_wrong_terminating_breakpoints = 0
	total_real_breakpoints = 0
	total_found_real_breakpoints = 0

	try:
		with open(args.compare_breakpoints, 'r') as cb_file:
			compare_breakpoints = json.load(cb_file)
	except Exception as e:
		logger.error("Could not read compare breakpoints file: '%s'", args.compare_breakpoints)
		print(e)
		sys.exit(1)

	for breakpoints in compare_breakpoints.values():
		breakpoints.sort()

	for contig_name in compare_breakpoints.keys():
		contig_breakpoints.setdefault(contig_name, [])

	for contig_name, breakpoints in contig_breakpoints.items():
		contig_len = len(contigs[contig_name].seq)
		last_existing_block = contig_len // args.blocking_size
		real_breakpoints = compare_breakpoints.setdefault(contig_name, [])
		matched_breakpoints = [False] * len(breakpoints)
		matched_real_breakpoints = [False] * len(real_breakpoints)

		for index, (breakpoint_start, breakpoint_end) in enumerate(breakpoints):
			start = breakpoint_start * args.blocking_size
			end = (breakpoint_end + 1) * args.blocking_size

			# Search overlapping real breakpoints
			for real_index, (real_start, real_end) in enumerate(real_breakpoints):
				if real_start <= end and start <= real_end:
					matched_breakpoints[index] = True
					matched_real_breakpoints[real_index] = True

		if len(breakpoints) > 0:
			if not matched_breakpoints[0] and breakpoints[0][0] == 0:
				total_wrong_terminating_breakpoints += 1
			if not matched_breakpoints[-1] and breakpoints[-1][1] == last_existing_block:
				total_wrong_terminating_breakpoints += 1

		total_correct_breakpoints += matched_breakpoints.count(True)
		total_wrong_breakpoints += matched_breakpoints.count(False)
		total_real_breakpoints += len(real_breakpoints)
		total_found_real_breakpoints += matched_real_breakpoints.count(True)

	logger.info("Total correct breakpoints: %d", total_correct_breakpoints)
	logger.info("Total wrong breakpoints: %d", total_wrong_breakpoints)
	logger.info("Total wrong non-terminating breakpoints: %d", total_wrong_breakpoints - total_wrong_terminating_breakpoints)
	logger.info("Total real breakpoints: %d", total_real_breakpoints)
	logger.info("Total missing real breakpoints: %d", total_real_breakpoints - total_found_real_breakpoints)
	logger.info("False positive rate: %.2f", total_wrong_breakpoints / total_breakpoints)
	logger.info("False positive rate without non-terminating: %.2f", (total_wrong_breakpoints - total_wrong_terminating_breakpoints) / total_breakpoints)
	logger.info("False negative rate: %.2f", 1.0 - (total_found_real_breakpoints / total_real_breakpoints))


### Split contigs ###

logger.info("Splitting contigs and writing them to: '%s'", args.output_contigs)
total_broken_contigs = 0

def generate_broken_contigs():
	global total_broken_contigs

	for contig_name, contig in contigs.items():
		breakpoints = contig_breakpoints[contig_name]

		if contig_name in unmapped_contigs:
			continue

		if len(breakpoints) == 0:
			total_broken_contigs += 1
			yield contig
			continue

		contig_start = 0
		contig_index = 0
		for (breakpoint_block_start, breakpoint_block_end) in breakpoints:
			breakpoint_start = args.blocking_size * breakpoint_block_start
			breakpoint_limit = args.blocking_size * (breakpoint_block_end + 1)

			if breakpoint_start - contig_start >= args.min_broken_contig_len:
				total_broken_contigs += 1
				yield SeqRecord(contig.seq[contig_start:breakpoint_start], id = contig_name + "." + str(contig_index))

			contig_start = breakpoint_limit
			contig_index += 1

SeqIO.write(generate_broken_contigs(), args.output_contigs, "fasta")

### Output final statistics ###

logger.info("Finished, printing statistics")
logger.info("Total input contigs: %d", len(contigs))
logger.info("Unmapped contigs: %d", len(unmapped_contigs))
logger.info("Unmapped reads: %d", unmapped_reads)
logger.info("Total breakpoints: %d", total_breakpoints)
logger.info("Total non-terminating breakpoints: %d", non_terminating_breakpoints)
logger.info("Total output contigs: %d", total_broken_contigs)
logger.info("Average evidence per block: %.2f", total_evidence / total_blocks)