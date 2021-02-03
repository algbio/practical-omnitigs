#!/usr/bin/env python
import subprocess
import sys

jobid = None
for arg in sys.argv[1:]:
	if arg.isdigit():
		jobid = arg

if arg is None:
	sys.exit("Not numeric argument given.")

states = [state.strip() for state in subprocess.check_output("sacct -j %s --format State --noheader -M all | awk '{print $1}'" % jobid, shell=True).decode(sys.stdout.encoding).strip().split('\n')]

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED", "REVOKED", "REQUEUED", "RESIZING"]
all_completed = True
for state in states:
	if state != "COMPLETED":
		all_completed = False

if all_completed:
  print("success")
elif any(r in states for r in running_status):
  print("running")
else:
  print("failed")