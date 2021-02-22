#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import os
import sys
import errno
import shutil

def silentremove(filename):
    try:
        if isfile(filename):
            os.remove(filename)
        else:
            shutil.rmtree(filename)
        print("Removed {}".format(filename))
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
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
		has_outfiles = False

		try:
			outfiles = alllines.split("Select jobs to execute...")[1]
			outfiles = outfiles.split("output:")[1]
			outfiles = outfiles.split("jobid:")[0].strip()
			outfiles = outfiles.split(", ")
			has_outfiles = True
		except IndexError:
			# Apparently, this job did not even get to output its metadata, so it surely produced no output files.
			pass

		if has_outfiles and "1 of 1 steps (100%) done" not in alllines:
			for out_file in outfiles:
				silentremove(out_file)

