#!/usr/bin/python3

"""
Compile a report for the HamCircuit pipeline.
Arguments: <file prefix>
"""

import sys

file_name_prefix = sys.argv[1]
maxn = int(sys.argv[2])
amount = maxn + 1

hamiltonian_graphs = 0
raw_tsp_runtime = 0.0
preprocessed_tsp_runtime = 0.0

for i in range(amount):
	file = open(file_name_prefix + "-" + str(i) + ".report", 'r')
	lines = file.readlines()
	file.close()

	for line in lines:
		if "Graph is hamiltonian" in line:
			hamiltonian_graphs += 1

		if "Raw TSP Runtime: " in line:
			raw_tsp_runtime += float(line.strip()[16:])

		if "Preprocessed TSP Runtime: " in line:
			preprocessed_tsp_runtime += float(line.strip()[25:])

report_file = open(file_name_prefix + ".0-" + str(maxn) + ".overallreport", 'w')

report_file.write("Hamiltonian fraction: " + str(float(hamiltonian_graphs) / float(amount) * 100.0) + "%\n")
report_file.write("Average raw runtime: " + str(raw_tsp_runtime / amount) + "\n")
report_file.write("Average preprocessed runtime: " + str(preprocessed_tsp_runtime / amount) + "\n")
