#!/usr/bin/python3

"""
Compile a report for the HamCircuit pipeline.
Arguments: <file prefix>
"""

import sys

file_name_prefix = sys.argv[1]
maxn = int(sys.argv[2])
n = int(sys.argv[3])
c = float(sys.argv[4])
amount = maxn + 1

hamiltonian_graphs = 0
raw_tsp_runtime = 0.0
preprocessed_tsp_runtime = 0.0
preprocessing_solved = 0
preprocessing_removed_nodes = 0.0
preprocessing_removed_edges = 0.0

for i in range(amount):
	file = open(file_name_prefix + "-" + str(i) + ".n" + str(n) + "-c" + str(c) + ".report", 'r')
	lines = file.readlines()
	file.close()

	for line in lines:
		if "Preprocessing the graph revealed that it is not hamiltonian" in line:
			preprocessing_solved += 1

		if "Preprocessing removed" in line:
			parts = line.strip().split("%")
			preprocessing_removed_nodes += float(parts[0].strip().split(" ")[-1].strip())
			preprocessing_removed_edges += float(parts[1].strip().split(" ")[-1].strip())

		if "Graph is hamiltonian" in line:
			hamiltonian_graphs += 1

		if "Raw TSP Runtime: " in line:
			raw_tsp_runtime += float(line.strip()[16:])

		if "Preprocessed TSP Runtime: " in line:
			preprocessed_tsp_runtime += float(line.strip()[25:])

report_file = open(file_name_prefix + ".0-" + str(maxn) + ".n" + str(n) + "-c" + str(c) + ".overallreport", 'w')

report_file.write("Hamiltonian fraction: " + str(float(hamiltonian_graphs) / float(amount) * 100.0) + "%\n")
report_file.write("Average raw runtime: " + str(raw_tsp_runtime / amount) + "\n")
report_file.write("Average preprocessed runtime: " + str(preprocessed_tsp_runtime / amount) + "\n")
report_file.write("Preprocessing solved cases: " + str(float(preprocessing_solved) / float(amount) * 100.0) + "%\n")
report_file.write("Preprocessing removed average nodes in unsolved cases: " + str(preprocessing_removed_nodes / (amount - preprocessing_solved)) + "%\n")
report_file.write("Preprocessing removed average edges in unsolved cases: " + str(preprocessing_removed_edges / (amount - preprocessing_solved)) + "%\n")