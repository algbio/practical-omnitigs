#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import os
import sys
import errno
import shutil

log_path = "logs/latest"

if len(sys.argv) >= 2:
	log_path = sys.argv[1]

try:
	log_files = [join(log_path, f) for f in listdir(log_path) if isfile(join(log_path, f))]
except Exception as e:
	print("Could not list log files")
	print(e)
	sys.exit(0)

for log_file in log_files:
	if "practical-omnitigs" not in log_file:
		continue

	with open(log_file, 'r') as f:
		alllines = '\n'.join(f.readlines())

		if "1 of 1 steps (100%) done" not in alllines:
			print(log_file)

