#!/usr/bin/python3

"""
Convert the output of the different validation tools into a LaTeX file.
Arguments: <ContigValidator results> <quast report.tex> <cli verify LaTeX results> <bandage png> <output file>
"""

import sys

contig_validator_file_name = sys.argv[1]
quast_file_name = sys.argv[2]
graph_statistics_file_name = sys.argv[3]
bandage_png_name = sys.argv[4]
output_file_name = sys.argv[5]

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

#########################################
### Process CLI graph statistics file ###
#########################################

graph_statistics_file = open(graph_statistics_file_name, 'r')
graph_statistics_lines = graph_statistics_file.readlines()
graph_statistics_lines = [x.strip() for x in graph_statistics_lines]

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

def write_image(output_file, caption, name, natwidth, natheight):
	pixel_pt_factor = 0.7
	output_file.write("\\begin{figure*}\n")
	output_file.write("\\centering\n")
	output_file.write("\\includegraphics[width=\\textwidth,natwidth=" + str(natwidth * pixel_pt_factor) + "pt,natheight=" + str(natheight * pixel_pt_factor) + "pt]{" + name + "}\n")
	output_file.write("\\end{figure*}\n")

output_file = open(output_file_name, 'w')
output_file.write(
	"""
	\\documentclass[10pt,a4paper,twocolumn]{article}
	\\usepackage[cm]{fullpage}
	\\usepackage{graphicx}
	\\begin{document}
	\\begin{description}
		\\item[Attention:] this file was produced automatically, and some statistics might not make sense for certain pipelines.
	\\end{description}
	"""
)

write_table(output_file, "Genome Graph Statistics", graph_statistics_lines)

write_table(output_file, "ContigValidator", contig_validator_lines)

write_table(output_file, "QUAST: \\# of contigs", quast_lines[0:6])
write_table(output_file, "QUAST: total length of contigs", quast_lines[6:12])
write_table(output_file, "QUAST: statistics for contigs $\\geq$ 500bp", quast_lines[12:26])
write_table(output_file, "QUAST: alignment statistics for contigs $\\geq$ 500bp", quast_lines[26:])

output_file.write("\\newpage")
write_image(output_file, "Bandage genome graph", bandage_png_name, 1000, 1000)


output_file.write(
	"""
	\\end{document}
	"""
)