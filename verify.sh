#!/bin/bash

in_file="GCF_000008865.2_ASM886v2_genomic.fna.unitigs.fa"
out_file="$in_file.verify"
clean_in_file="$in_file.clean"

cargo run -- --input "$in_file" --output "$out_file" verify

function join_by { local IFS="$1"; shift; echo "$*"; }
while read -r line
do
	if [[ $line == \>* ]]; then
		arguments=($line)
		result="${arguments[0]} ${arguments[1]} ${arguments[2]} ${arguments[3]}"
		arguments=("${arguments[@]:4}")
		IFS=$'\n' arguments=($(sort -t ":" -k2,2r -k3,3n <<<"${arguments[*]}"))
		unset IFS
		arguments=$(printf " %s" "${arguments[@]}")
		result=`echo "$result$arguments" | sed 's/ *$//g'`
		echo "${result}"
	else
		echo "$line"
	fi
done < "$in_file" > "$clean_in_file"


diff=`diff "$clean_in_file" "$out_file"`
if [ -n "$diff" ]; then
	echo "$diff"
	exit 1
fi
