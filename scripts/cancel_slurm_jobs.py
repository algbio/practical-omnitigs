#!/usr/bin/env python
import subprocess
import sys
import traceback

jobids = []
for arg in sys.argv[1:]:
	if arg.isdigit():
		jobids.append(arg)

if len(jobids) == 0:
	sys.exit("No numeric argument given.")

output = subprocess.check_output(f"scancel {" ".join(jobids)}", shell = True).decode(sys.stdout.encoding).strip()
print(output)