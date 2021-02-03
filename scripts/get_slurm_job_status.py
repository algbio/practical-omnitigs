#!/usr/bin/env python
import subprocess
import sys

jobid = None
for arg in sys.argv[1:]:
	if arg.isdigit():
		jobid = arg

if arg is None:
	sys.exit("Not numeric argument given.")

output = str(subprocess.check_output("sacct -j %s --format State --noheader -M all | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED", "REVOKED", "REQUEUED", "RESIZING"]
if "COMPLETED" in output:
  print("success")
elif any(r in output for r in running_status):
  print("running")
else:
  print("failed")