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


contig_validator_file = open(contig_validator_file_name, 'r')
contig_validator_lines = contig_validator_file.readlines()
contig_validator_lines = contig_validator_lines[1].split()[1:] # Remove header and file name column
contig_validator_lines[0] = "\\%exact & " + contig_validator_lines[0] + "\\%\\\\"
contig_validator_lines[1] = "\\%align & " + contig_validator_lines[1].replace("%", "\\%") + "\\\\"
contig_validator_lines[2] = "recall & " + contig_validator_lines[2].replace("%", "\\%") + "\\\\"
contig_validator_lines[3] = "precision & " + contig_validator_lines[3].replace("%", "\\%") + "\\\\"


########################
### Build LaTeX file ###
########################

def table_header(caption):
	return """
	\\begin{table}[ht]
	\\begin{center}
	\\caption{""" + caption + """}
	\\begin{tabular}{|l*{1}{|r}|}
	\\hline
	Parameter & Value \\\\ \\hline
	"""

table_footer = """\\hline
	\\end{tabular}
	\\end{center}
	\\end{table}
	"""

def write_table(output_file, caption, rows):
	output_file.write(table_header(caption))
	for row in rows:
		output_file.write(row + '\n')
	output_file.write(table_footer)

output_file = open(output_file_name, 'w')
output_file.write(
	"""
	\\documentclass[10pt,a4paper,twocolumn]{article}
	\\usepackage[cm]{fullpage}
	\\begin{document}
	\\begin{description}
		\\item[Attention:] this file was produced automatically, and some statistics might not make sense for certain pipelines.
	\\end{description}
	"""
)
#output_file.write(table_header("QUAST: All statistics are based on contigs of size $\\geq$ 500 bp, unless otherwise noted (e.g., \"\\# contigs ($\\geq$ 0 bp)" and "Total length ($\\geq$ 0 bp)\" include all contigs)\\newline	ContigValidator: Minor errors due to circularisation may occur."))
write_table(output_file, "QUAST: \\# of contigs", quast_lines[0:6])
write_table(output_file, "QUAST: total length of contigs", quast_lines[6:12])
write_table(output_file, "QUAST: statistics for contigs $\\geq$ 500bp", quast_lines[12:26])
write_table(output_file, "QUAST: alignment statistics for contigs $\\geq$ 500bp", quast_lines[26:])


write_table(output_file, "ContigValidator", contig_validator_lines)

output_file.write(
	"""
	\\end{document}
	"""
)