#!/usr/bin/python3

"""
Convert the output of the different validation tools into a LaTeX file.
Arguments: <genome name> <cli verify LaTeX results> <bandage png> <output file> [<label> <experiment prefix> ...] 
"""

import sys, subprocess, hashlib, os, pathlib
import traceback

hashdir = os.path.abspath(sys.argv[1])
genome_name_file_name = sys.argv[2]
graph_statistics_file_name = sys.argv[3]
bandage_png_name = sys.argv[4]
combined_eaxmax_plot_name = sys.argv[5]
output_file_name = sys.argv[6]

experiments = []

for i in range(7, len(sys.argv), 2):
	if i + 1 >= len(sys.argv):
		print("Uneven number of experiment parameters. Each parameter must be a pair of <label> and <experiment prefix>")
		exit(1)

	experiments.append((sys.argv[i], sys.argv[i + 1]))

def append_latex_table_second_column(table, appendix):
	if len(table) == 0:
		return appendix

	appendix_map = {}
	for line in appendix:
		key = line[:line.index("&")].strip()
		value = line[line.index("&") - 1:]
		appendix_map[key] = value

	result = []
	for line in table:
		line = line.strip()
		if line[-2:] == "\\\\":
			line = line[:-2] # Remove trailing new line backslashes

		key = line[:line.index("&")].strip()
		appendix_value = appendix_map.pop(key, None)

		if appendix_value is None:
			print("Warning: Missing value {} in appendix".format(key))
			result.append(line + " & N/A \\\\")
		else:
			result.append(line + appendix_value)

	if len(appendix_map) != 0:
		sys.exit("Some appended metrices were not in the table before: {}".format(appendix_map))

	return result

#########################
### Process name file ###
#########################

name_file = open(genome_name_file_name, 'r')
name_lines = name_file.readlines()
name_lines = [x.replace("_", "\\_") for x in name_lines]

###############################
### Process algorithm files ###
###############################

def read_algorithm_file(prefix):
	try:
		algorithm_file = open(os.path.join(prefix, "compute_injectable_contigs.tex"))
		algorithm_lines = algorithm_file.readlines()
		return algorithm_lines
	except:
		print("Did not find algorithm file for '" + prefix + "'")
		return []

headline = "Parameter"
algorithm_table = []
for label, prefix in experiments:
	headline += " & " + label
	table = read_algorithm_file(os.path.join(prefix, "injectable_contigs"))
	algorithm_table = append_latex_table_second_column(algorithm_table, table)

algorithm_table = [headline + "\\\\ \\hline"] + algorithm_table


##########################
### Process QUAST files ###
##########################

def read_quast_file(filename):
	try:
		quast_file = open(filename, 'r')
		quast_lines = quast_file.readlines()
		quast_lines = [x.strip() for x in quast_lines]
		quast_lines = [x.replace("\\hline", "") for x in quast_lines]
		quast_lines = quast_lines[8:-4] # Remove LaTeX header and footer
		quast_lines = [line.replace("misassemblies caused by fragmented reference", "mis. caused by frag. ref.") for line in quast_lines]
		return quast_lines
	except:
		print("Did not find quast file for '" + filename + "'")
		return []

headline = "Parameter"
quast_table = []
for (label, prefix) in experiments:
	headline += " & " + label
	table = read_quast_file(prefix + "quast/report.tex")
	quast_table = append_latex_table_second_column(quast_table, table)

quast_table = [headline + "\\\\ \\hline"] + quast_table

headline = "Parameter"
quast_misassemblies_table = []
for (label, prefix) in experiments:
	headline += " & " + label
	table = read_quast_file(prefix + "quast/contigs_reports/misassemblies_report.tex")
	quast_misassemblies_table = append_latex_table_second_column(quast_misassemblies_table, table)

quast_misassemblies_table = [headline + "\\\\ \\hline"] + quast_misassemblies_table


####################################
### Process ContigValidator files ###
####################################

def read_contig_validator_file(prefix):
	try:
		contig_validator_file = open(prefix + ".contigvalidator", 'r')
		contig_validator_lines = contig_validator_file.readlines()
		contig_validator_lines = contig_validator_lines[1].split()[1:] # Remove header and file name column
		contig_validator_lines[0] = "\\%exact & " + contig_validator_lines[0] + "\\%\\\\"
		contig_validator_lines[1] = "\\%align & " + contig_validator_lines[1].replace("%", "\\%") + "\\\\"
		contig_validator_lines[2] = "recall & " + contig_validator_lines[2].replace("%", "\\%") + "\\\\"
		contig_validator_lines[3] = "precision & " + contig_validator_lines[3].replace("%", "\\%") + "\\\\"
		return contig_validator_lines
	except:
		print("Did not find contig validator file for '" + prefix + "'")
		return []

headline = "Parameter"
contig_validator_table = []
for (label, prefix) in experiments:
	headline += " & " + label
	table = read_contig_validator_file(prefix)
	contig_validator_table = append_latex_table_second_column(contig_validator_table, table)

contig_validator_table = [headline + "\\\\ \\hline"] + contig_validator_table

##########################################
### Process CLI graph statistics files ###
##########################################

try:
	graph_statistics_file = open(graph_statistics_file_name, 'r')
	graph_statistics_lines = graph_statistics_file.readlines()
	graph_statistics_table = ["Parameter & Value \\\\ \\hline"] + [x.strip() for x in graph_statistics_lines]
except:
	print("Did not find graph statistics file '" + graph_statistics_file_name + "'")
	graph_statistics_table = []

########################
### Build LaTeX file ###
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
	\\fontsize{6pt}{7pt}\\selectfont
	\\caption{\\fontsize{6pt}{7pt}\\selectfont\\ """ + caption + """}
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

def write_image(output_file, caption, file, natwidth, natheight):
	hasher = hashlib.sha3_512()
	hasher.update(file.encode())
	hashlinkdir = os.path.join(hashdir, str(hasher.hexdigest()))
	hashlink = os.path.join(hashlinkdir, "file." + file.split(".")[-1])

	try:
		os.remove(hashlink)
	except OSError:
		pass

	pathlib.Path(hashlinkdir).mkdir(parents=True, exist_ok=True)

	try:
		os.symlink(os.path.abspath(file), hashlink)
	except FileExistsError:
		print("error: could not create symlink because it already exists. We assume that it is correct.")
	except OSError:
		print("Error: could not create symlink, but just continuing because it might have been correctly created concurrently.")
		traceback.print_exc()

	pixel_pt_factor = 0.7
	output_file.write("\\begin{figure*}\n")
	output_file.write("\\centering\n")
	output_file.write("\\includegraphics[width=\\textwidth,natwidth=" + str(natwidth * pixel_pt_factor) + "pt,natheight=" + str(natheight * pixel_pt_factor) + "pt]{" + str(hashlink) + "}\n")
	output_file.write("\\caption{" + str(caption) + "}")
	output_file.write("\\end{figure*}\n")

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

for line in name_lines:
	output_file.write("\\item " + line)
output_file.write("\\end{itemize}\n")

write_table(output_file, "Genome Graph Statistics", 1, graph_statistics_table)

write_table(output_file, "Algorithm Statistics", len(experiments), algorithm_table)

write_table(output_file, "ContigValidator", len(experiments), contig_validator_table)

write_table(output_file, "QUAST: \\# of contigs", len(experiments), quast_table[0:7])
write_table(output_file, "QUAST: total length of contigs", len(experiments), [quast_table[0]] + quast_table[7:13])
write_table(output_file, "QUAST: statistics for contigs $\\geq$ 500bp (or 3000bp for QUAST-LG)", len(experiments), [quast_table[0]] + quast_table[13:28])
write_table(output_file, "QUAST: alignment statistics for contigs $\\geq$ 500bp (or 3000bp for QUAST-LG)", len(experiments), [quast_table[0]] + quast_table[28:])
write_table(output_file, "QUAST: misassembly statistics for contigs $\\geq$ 500bp (or 3000bp for QUAST-LG)", len(experiments), quast_misassemblies_table)


from os import path

if combined_eaxmax_plot_name != 'none':
    output_file.write("\\newpage")
    write_image(output_file, "EAxmax", combined_eaxmax_plot_name, 1000, 1000)

output_file.write("\\newpage")
for label, prefix in experiments:
	plot_file_path = prefix + "quast/aligned_stats/EAxmax_plot.pdf"
	if path.isfile(plot_file_path):
		quast_png_name = plot_file_path
		write_image(output_file, "QUAST EAxmax graph for " + label, quast_png_name, 1000, 1000)
	else:
		print("Did not find '" + plot_file_path + "'")

for label, prefix in experiments:
	plot_file_path = prefix + "quast/aligned_stats/NGAx_plot.pdf"
	if path.isfile(plot_file_path):
		quast_png_name = plot_file_path
		write_image(output_file, "QUAST NGAx graph for " + label, quast_png_name, 1000, 1000)
	else:
		print("Did not find '" + plot_file_path + "'")

if bandage_png_name != 'none':
	output_file.write("\\newpage")
	write_image(output_file, "Bandage genome graph", bandage_png_name, 1000, 1000)


output_file.write(
	"""
	\\end{landscape}
	\\end{document}
	"""
)