#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import errno

def silentremove(filename):
    try:
        os.remove(filename)
        print("Removed {}".format(filename))
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
log_path = "logs/latest"
log_files = [join(log_path, f) for f in listdir(log_path) if isfile(join(log_path, f))]

for log_file in log_files:
	if "practical-omnitigs" not in log_file:
		continue

	with open(log_file, 'r') as f:
		alllines = '\n'.join(f.readlines())

		outfiles = alllines.split("Select jobs to execute...")[1]
		outfiles = outfiles.split("output:")[1]
		outfiles = outfiles.split("jobid:")[0].strip()
		outfiles = outfiles.split(", ")

		if "1 of 1 steps (100%) done" not in alllines:
			for out_file in outfiles:
				silentremove(out_file)

