#!/usr/bin/python3

"""
Convert the output of the different validation tools into a LaTeX file.
Arguments: <ContigValidator results> <quast report.tex> <output file>
"""

import sys

contig_validator_file_name = sys.argv[1]
quast_file_name = sys.argv[2]
output_file_name = sys.argv[3]

##########################
### Process QUAST file ###
##########################

quast_file = open(quast_file_name, 'r')
quast_lines = quast_file.readlines()
quast_lines = [x.strip() for x in quast_lines]
quast_lines = [x.replace("\\hline", "") for x in quast_lines]
quast_lines = quast_lines[8:-4] # Remove LaTeX header and footer

####################################
### Process ContigValidator file ###
####################################

# TODO
contig_validator_lines = []

########################
### Build LaTeX file ###
########################

output_file = open(output_file_name, 'w')
output_file.write(
	"""
	\\documentclass[10pt,a4paper]{article}
	\\usepackage{fullpage}
	\\begin{document}
	\\begin{table}[ht]
	\\begin{center}
	\\caption{QUAST: All statistics are based on contigs of size $\\geq$ 500 bp, unless otherwise noted (e.g., "\\# contigs ($\\geq$ 0 bp)" and "Total length ($\\geq$ 0 bp)" include all contigs)\\newline
	ContigValidator: Minor errors due to circularisation may occur.}
	\\begin{tabular}{|l*{1}{|r}|}
	\\hline
	Parameter & Value \\\\
	"""
)

output_file.write("\\textbf{QUAST} & \\\\")
for line in quast_lines:
	output_file.write(line + '\n')


output_file.write("\\textbf{ContigValidator} & \\\\")
for line in contig_validator_lines:
	output_file.write(line + '\n')

output_file.write(
	"""
	\\end{tabular}
	\\end{center}
	\\end{table}
	\\end{document}
	"""
)