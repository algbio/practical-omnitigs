#!/usr/bin/env python
import subprocess
import sys

jobid = None
for arg in sys.argv[1:]:
	if arg.isdigit():
		jobid = arg

if arg is None:
	sys.exit("Not numeric argument given.")

states = [state.strip() for state in subprocess.check_output("sacct -j %s --format 'JobID%20,State%20' --noheader -M all" % jobid, shell=True).decode(sys.stdout.encoding).strip().split('\n')]
states = [(line.split()[0], line.split()[1]) for line in states]
states = [state for id, state in states if id.isdigit()]


running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED", "REVOKED", "REQUEUED", "RESIZING"]
completed_status = ["COMPLETED"]

if any(r in states for r in completed_status):
  print("success")
elif any(r in states for r in running_status):
  print("running")
else:
  print("failed")