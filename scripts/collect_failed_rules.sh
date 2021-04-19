#!/bin/bash

set -e

if [ -z "$1" ]; then
	echo "Error: require the log directory as first argument"
	exit 1
fi

JOBS_LOG="logs/$1/jobs.log"
SEFF_LOG="logs/$1/jobs_seff.log"
FAILED_JOBS_LOG="logs/$1/jobs_failed_seff.log"

if [[ ! -f "$JOBS_LOG" ]]; then
	echo "Error: input file does not exist: $JOBS_LOG"
	exit 1
fi

rm -f "$SEFF_LOG" "$FAILED_JOBS_LOG"
touch "$SEFF_LOG" "$FAILED_JOBS_LOG"

JOB_IDS=$(grep -oE "^JobId=[0-9]{3,12}" "$JOBS_LOG" | sed "s/JobId=//g")
JOBS_HAVE_FAILED=""
JOBS_NOT_FOUND=""

set +e

for JOB_ID in $JOB_IDS; do
	if [ ! -z "$(seff $JOB_ID 2>&1 | grep 'Job not found.')" ]; then
		echo "Job ID: $JOB_ID" >> $FAILED_JOBS_LOG
		echo "Job ID: $JOB_ID" >> $SEFF_LOG
		seff $JOB_ID >> $FAILED_JOBS_LOG 2>&1
		echo "" >> $FAILED_JOBS_LOG
		echo "Job $JOB_ID cannot be found anymore."
		JOBS_NOT_FOUND="true"
	elif [ -z "$(seff $JOB_ID 2>&1 | grep '(exit code 0)')" ]; then
		seff $JOB_ID >> $FAILED_JOBS_LOG 2>&1
		echo "" >> $FAILED_JOBS_LOG
		echo "Job $JOB_ID has failed."
		JOBS_HAVE_FAILED="true"
	fi

	seff $JOB_ID >> $SEFF_LOG 2>&1
	echo "" >> $SEFF_LOG
done

echo "Wrote all seffs to $SEFF_LOG."

if [ "$JOBS_HAVE_FAILED" == "true" ]; then
	echo "Some jobs have failed. See $FAILED_JOBS_LOG for more details."
fi

if [ "$JOBS_NOT_FOUND" == "true" ]; then
	echo "Some jobs could not be found. See $FAILED_JOBS_LOG for more details."
fi
