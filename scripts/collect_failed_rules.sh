#!/bin/bash

set -e

if [ -z "$1" ]; then
	echo "Error: require the log directory as first argument"
	exit 1
fi

JOBS_LOG="logs/$1/jobs.log"
SEFF_LOG="logs/$1/jobs_seff.log"
FAILED_JOBS_LOG="logs/$1/failed_jobs_seff.log"

if [[ ! -f "$JOBS_LOG" ]]; then
	echo "Error: input file does not exist: $JOBS_LOG"
	exit 1
fi

JOB_IDS=$(grep -oE "^JobId=[0-9]{3,12}" "$JOBS_LOG" | sed "s/JobId=//g")
JOBS_HAVE_FAILED=""

for JOB_ID in $JOB_IDS; do
	if [ -z "$(seff $JOB_ID | grep '(exit code 0)')" ]; then
		seff $JOB_ID >> $FAILED_JOBS_LOG
		echo "" >> $FAILED_JOBS_LOG
		echo "Job $JOB_ID has failed."
		JOBS_HAVE_FAILED="true"
	fi

	seff $JOB_ID >> $SEFF_LOG
	echo "" >> $SEFF_LOG
done

echo "Wrote all seffs to $SEFF_LOG."

if [ "$JOBS_HAVE_FAILED" == "true" ]; then
	echo "Jobs have failed. See $FAILED_JOBS_LOG for more details."
fi
