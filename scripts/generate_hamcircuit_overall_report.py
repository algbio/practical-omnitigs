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
preprocessing_runtime = 0.0
preprocessed_tsp_runtime = 0.0
preprocessed_total_runtime = 0.0
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

		if "Preprocessing Runtime: " in line:
			preprocessing_runtime += float(line.strip()[22:])

		if "Preprocessed TSP Runtime: " in line:
			preprocessed_tsp_runtime += float(line.strip()[25:])

		if "Preprocessed Total Runtime: " in line:
			preprocessed_total_runtime += float(line.strip()[27:])

report_file = open(file_name_prefix + ".0-" + str(maxn) + ".n" + str(n) + "-c" + str(c) + ".overallreport", 'w')

report_file.write("Hamiltonian fraction: {:.1f}%\n".format(float(hamiltonian_graphs) / float(amount) * 100.0))
report_file.write("Average preprocessing runtime: {:.1f} seconds\n".format(preprocessing_runtime / amount))
report_file.write("Average preprocessed TSP runtime: {:.1f} seconds\n".format(preprocessed_tsp_runtime / amount))
report_file.write("Average preprocessed total runtime: {:.1f} seconds\n".format(preprocessed_total_runtime / amount))
report_file.write("Average raw TSP runtime: {:.1f} seconds\n".format(raw_tsp_runtime / amount))
report_file.write("Preprocessing solved cases: {:.1f}%\n".format(float(preprocessing_solved) / float(amount) * 100.0))
report_file.write("Preprocessing removed average nodes in unsolved cases: {:.1f}%\n".format(preprocessing_removed_nodes / (amount - preprocessing_solved)))
report_file.write("Preprocessing removed average edges in unsolved cases: {:.1f}%\n".format(preprocessing_removed_edges / (amount - preprocessing_solved)))