#!/usr/bin/env python3

import sys

nodes_path = sys.argv[1]

with open(nodes_path, 'r') as nodes_file:
	for line in nodes_file:
		columns = [s for s in line.split(" ") if len(s) > 0]

		aligns = []
		for column in columns[2:]:
			column = column.strip()
			print(column)
			fields = column.split("_")
			fields = ["_".join(fields[:-5])] + fields[-5:]

#NC_002695.2_3080719_3105661_F_41_4 -> 3091215_3092239
#NC_002695.2_3110074_3085294_R_69_4 -> 3091386_3092410 (shift 171)
#NC_002695.2_3115177_3090418_R_89_4 -> 3091369_3092393 (shift 154)
#NC_002695.2_3086190_3110151_F_20_4 -> 3091310_3092334 (shift 95)

			read_start = int(fields[1])
			read_end = int(fields[2])
			forward = fields[3] == "F"
			align_read_offset = int(fields[4]) * 256
			align_len = int(fields[5]) * 256
			if read_start < read_end:
				align_start = read_start + align_read_offset
				align_end = align_start + align_len
			else:
				read_start, read_end = read_end, read_start
				forward = not forward
				align_end = read_end - align_read_offset
				align_start = align_end - align_len


			aligns.append((fields[0], read_start, read_end, align_start, align_end, forward))

		print(aligns)