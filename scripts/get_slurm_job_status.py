#!/usr/bin/env python
import subprocess
import sys
import traceback

RUNNING_STATUS = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED", "REVOKED", "REQUEUED", "RESIZING"]
COMPLETED_STATUS = ["COMPLETED"]
FAILED_STATUS = ["FAILED", "CANCELLED", "BOOT_FAIL", "DEADLINE", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "TIMEOUT"]
ALL_STATUSES = RUNNING_STATUS + COMPLETED_STATUS + FAILED_STATUS

jobid = None
joined_args = "' '".join(sys.argv)
for arg in sys.argv[1:]:
  for subarg in arg.split():
    if subarg.isdigit():
      if jobid is not None:
        sys.exit(f"Found two numeric arguments.argument Arguments: '{joined_args}'")
      jobid = subarg

if jobid is None:
  sys.exit(f"No numeric argument given. Arguments: '{joined_args}'")

try:
  original_states = [state.strip() for state in subprocess.check_output("sacct -j {} --format 'JobID%20,State%20' --noheader -M all".format(jobid), shell=True).decode(sys.stdout.encoding).strip().split('\n')]
except Exception:
  traceback.print_exc(file = sys.stderr)
  print(f"Error running sacct\n", file = sys.stderr, flush = True)
  print("running")
  sys.exit(0)

if len(original_states) == 0:
  print(f"Error sacct\n{original_states}", file = sys.stderr, flush = True)
  print("running")
  sys.exit(0)

states = []

for line in original_states:
  state = None
  add = False
  for word in line.split():
    if word.isdigit():
      add = True
    if word in ALL_STATUSES:
      state = word

  if add:
    assert state is not None, f"Got a job id but no state in line: {line}"
    states.append(state)

if len(states) == 0:
  print(f"Error sacct\n{original_states}", file = sys.stderr, flush = True)
  print("running")
  sys.exit(0)

# TODO check with multi cluster submit
if all(r in COMPLETED_STATUS for r in states):
  print("success")
elif any(r in RUNNING_STATUS for r in states):
  print("running")
elif all(r in FAILED_STATUS for r in states):
  print("failed")
else:
  print("Unknown state combination: {}".format(states))
  sys.exit(1)