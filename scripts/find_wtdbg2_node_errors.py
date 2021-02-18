#!/usr/bin/env python3

import sys

nodes_path = sys.argv[1]

with open(nodes_path, 'r') as nodes_file:
	for line in nodes_file:
		columns = line.split("\t")
		print(columns[0])