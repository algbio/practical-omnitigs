#!/usr/bin/python3

import argparse
import sys

parser = argparse.ArgumentParser(description="Create an aggregated report of wtdbg2 experiments")
parser.add_argument("--source-reports", type=str, nargs='+', required=True)
parser.add_argument("--source-report-names", type=str, nargs='+', required=True)
parser.add_argument("--output", type=str, required=True)

args = parser.parse_args()

source_reports = args.source_reports
source_report_names = args.source_report_names
output_file_name = args.output

if len(source_reports) != len(source_report_names):
	sys.exit("Source reports and source report names differ in lengths")

########################
### Helper functions ###
########################

def table_header(caption, column_count):
	#header = """
	#\\begin{table}[ht]
	#\\begin{center}
	#\\caption{""" + caption + """}
	#\\begin{tabular}{|l*{1}{|r}|}
	#\\hline
	#"""
	header = """
	\\begin{table}[ht]
	\\begin{center}
	\\caption{""" + caption + """}
	\\begin{tabular}{|l*{1}|"""
	for _ in range(column_count):
		header += "r"
	header += """|}
	\\hline
	"""
	return header

table_footer = """\\hline
	\\end{tabular}
	\\end{center}
	\\end{table}
	"""

def write_table(output_file, caption, column_count, rows):
	output_file.write(table_header(caption, column_count))
	for row in rows:
		output_file.write(row + '\n')
	output_file.write(table_footer)

def write_aggregated_table(output_file, caption, metrics, metric_shortname=None):
	if type(metrics) is not list:
		metrics = [metrics]

	if metric_shortname is None:
		metric_shortname = '/'.join(metrics)

	global_value_lines = []
	global_table_header = None
	for source_report_name, source_report in zip(source_report_names, source_reports):
		with open(source_report, 'r') as input_file:
			source_report_lines = input_file.readlines()

		table_headers = []
		value_lines = []
		for metric in metrics:
			# Find applicable table
			last_table_header = None
			value_line = None
			for line in source_report_lines:
				if "Parameter" in line:
					last_table_header = line
				if metric in line:
					value_line = line

			if last_table_header is None or value_line is None:
				sys.exit("Found value line outside of table")
			table_headers.append(last_table_header)
			value_lines.append(value_line)

		# Combine metrics
		combined_value_line = None
		for table_header, value_line in zip(table_headers, value_lines):
			if global_table_header is None:
				global_table_header = table_header
			elif global_table_header != table_header:
				sys.exit("Found differing table headers: {} != {}".format(global_table_header, table_header))
			value_line = [value.strip() for value in value_line.split("&")[1:]]
			if combined_value_line is None:
				combined_value_line = value_line
			else:
				for index, value in enumerate(value_line):
					combined_value_line[index] += "/" + value
		global_value_lines.append(source_report_name.replace("_", "\\_") + " & " + " & ".join(combined_value_line) + "\\\\")

	# Add header and write table
	table = [global_table_header.replace("Parameter", metric_shortname.replace("_", "\\_"))] + global_value_lines
	write_table(output_file, caption, len(global_table_header.split("&")) - 1, table)

##############
### Header ###
##############

import subprocess
revision = subprocess.check_output(["git", "describe"]).strip()
output_file = open(output_file_name, 'w')

output_file.write(
	"""
	\\documentclass[12pt,a4paper]{article}
	\\usepackage[margin=0pt]{geometry}
	\\usepackage{lmodern}
	\\usepackage{pdflscape}
	\\usepackage[T1]{fontenc}
	\\usepackage{graphicx}
	\\begin{document}
	\\begin{landscape}
	\\fontsize{6pt}{7pt}\\selectfont
	\\begin{description}
		\\item[Attention:] this file was produced automatically, and some statistics might not make sense for certain pipelines.
		\\item[Revision:] """ + str(revision) + """
	\\end{description}
	This file contains statistics about the following genome(s):
	\\begin{itemize}
	"""
)

for source_report in source_reports:
	output_file.write("\\item " + source_report.replace("_", "\\_") + '\n')
output_file.write("\\end{itemize}\n")

#####################
### EAxmax tables ###
#####################

write_aggregated_table(output_file, "EA50max", "EA50max")

##########################
### Misassembly tables ###
##########################

write_aggregated_table(output_file, "Number of unique/total extensive misassemblies. These do not include local misassemblies. The `uniqueness' is determined heuristically.", ["# unique misassemblies", "# misassemblies"], "unique/total ext. mis.")
write_aggregated_table(output_file, "Number of unique/total local misassemblies. The `uniqueness' is determined heuristically.", ["# unique local misassemblies", "# local misassemblies"], "unique/total loc. mis.")

##############
### Footer ###
##############

output_file.write("\\end{landscape}\\end{document}")
