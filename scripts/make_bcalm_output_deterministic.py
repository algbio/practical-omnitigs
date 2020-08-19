#in_file="$1"
#clean_in_file="$2"

#while read -r line
#do
#	if [[ $line == \>* ]]; then
#		arguments=($line)
#		result="${arguments[0]} ${arguments[1]} ${arguments[2]} ${arguments[3]}"
#		arguments=("${arguments[@]:4}")
#		IFS=$'\n' arguments=($(sort -t ":" -k2,2r -k3,3n <<<"${arguments[*]}"))
#		unset IFS
#		arguments=$(printf " %s" "${arguments[@]}")
#		result=`echo "$result$arguments" | sed 's/ *$//g'`
#		echo "${result}"
#	else
#		echo "$line"
#	fi
#done < "$in_file" > "$clean_in_file"

import sys

input = open(sys.argv[1], 'r')
output = open(sys.argv[2], 'w')

def sort_key(token):
	bits = token.split(':')
	bits[1] = 'a' if bits[1] == '+' else 'b'
	bits[3] = 'a' if bits[3] == '-' else 'b'
	return (bits[1], int(bits[2]), bits[3])

for line in input:
	if line.startswith(">"):
		tokens = list(filter(None, line.split()))
		const_tokens = tokens[:4]
		tokens = tokens[4:]
		tokens = sorted(tokens, key = sort_key)
		line = ' '.join(const_tokens + tokens) + '\n'

	output.write(line)