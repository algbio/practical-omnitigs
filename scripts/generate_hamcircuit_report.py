#!/usr/bin/python3

"""
Compile a report for the HamCircuit pipeline.
Arguments: <file prefix>
"""

import sys

file_name_prefix = sys.argv[1]

preprocesslog_lines = open(file_name_prefix + ".preprocesslog", 'r').readlines()
solution_raw_lines = open(file_name_prefix + ".raw.sol", 'r').readlines()
solution_preprocessed_lines = open(file_name_prefix + ".preprocessed.sol", 'r').readlines()
tsplog_raw_lines = open(file_name_prefix + ".raw.tsplog", 'r').readlines()
tsplog_preprocessed_lines = open(file_name_prefix + ".preprocessed.tsplog", 'r').readlines()

report_file = open(file_name_prefix + ".report", 'w')

preprocessing_running_time = -1
for line in preprocesslog_lines:
	if "Preprocessing removed" in line:
		line = "Preprocessing removed" + line.strip().split("Preprocessing removed")[1]
		report_file.write(line + "\n")

	if "Preprocessing the graph revealed that it is not hamiltonian" in line:
		report_file.write("Preprocessing the graph revealed that it is not hamiltonian\n")

	if "Preprocessing took" in line:
		preprocessing_running_time = float(line.split("Preprocessing took")[1].split("seconds")[0].strip())

	if "Found" in line and "safe walks" in line and "safe walks of length > 2" not in line:
		report_file.write("Safe walks: " + line.strip().split("Found")[1].split("safe")[0].strip())

	if "Found" in line and "safe walks of length > 2" in line:
		report_file.write("Safe walks length > 2: " + line.strip().split("Found")[1].split("safe")[0].strip())

	if "Found" in line and "non-trivial safe walks" in line:
		report_file.write("Non-trivial safe walks: " + line.strip().split("Found")[1].split("non-trivial")[0].strip())

raw_number_of_nodes = -1
raw_optimal_solution = -1
raw_tsp_running_time = -1
for line in tsplog_raw_lines:
	if "Number of Nodes: " in line:
		line = line.strip()[16:]
		raw_number_of_nodes = int(round(float(line)))

	if "Optimal Solution: " in line:
		line = line.strip()[17:]
		raw_optimal_solution = int(round(float(line)))

	if "Total Running Time: " in line and "(seconds)" in line:
		line = line.strip()[19:].strip().split(" ")[0].strip()
		raw_tsp_running_time = float(line)

preprocessed_number_of_nodes = -1
preprocessed_optimal_solution = -1
preprocessed_tsp_running_time = -1
for line in tsplog_preprocessed_lines:
	if "Number of Nodes: " in line:
		line = line.strip()[16:]
		preprocessed_number_of_nodes = int(round(float(line)))

	if "Optimal Solution: " in line:
		line = line.strip()[17:]
		preprocessed_optimal_solution = int(round(float(line)))

	if "Total Running Time: " in line and "(seconds)" in line:
		line = line.strip()[19:].strip().split(" ")[0].strip()
		preprocessed_tsp_running_time = float(line)

if preprocessing_running_time == -1:
	print("Missing number preprocessing_running_time")
	sys.exit(1)

if raw_number_of_nodes == -1:
	print("Missing number raw_number_of_nodes")
	sys.exit(1)

if raw_optimal_solution == -1:
	print("Missing number raw_optimal_solution")
	sys.exit(1)

if raw_tsp_running_time == -1:
	print("Missing number raw_tsp_running_time")
	sys.exit(1)

if preprocessed_number_of_nodes == -1:
	print("Missing number preprocessed_number_of_nodes")
	sys.exit(1)

if preprocessed_optimal_solution == -1:
	print("Missing number preprocessed_optimal_solution")
	sys.exit(1)

if preprocessed_tsp_running_time == -1:
	print("Missing number preprocessed_tsp_running_time")
	sys.exit(1)

raw_hamiltonian = raw_number_of_nodes * 5 == raw_optimal_solution
preprocessed_hamiltonian = preprocessed_number_of_nodes * 5 == preprocessed_optimal_solution

if raw_hamiltonian != preprocessed_hamiltonian:
	print("Difference between hamiltonianess of raw and preprocessed graph.")
	report_file.write("Error: Difference in hamiltonianess\n")
	sys.exit(1)
else:
	report_file.write("Hamiltonianess matches\n")

if raw_hamiltonian:
	report_file.write("Graph is hamiltonian\n")
else:
	report_file.write("Graph is not hamiltonian\n")

report_file.write("Raw TSP Runtime: " + str(raw_tsp_running_time) + "\n")
report_file.write("Preprocessed TSP Runtime: " + str(preprocessed_tsp_running_time) + "\n")
report_file.write("Preprocessing Runtime: " + str(preprocessing_running_time) + "\n")
report_file.write("Preprocessed Total Runtime: " + str(preprocessing_running_time + preprocessed_tsp_running_time) + "\n")
