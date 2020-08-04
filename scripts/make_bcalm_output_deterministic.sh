#!/bin/bash

in_file="$1"
clean_in_file="$2"

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

