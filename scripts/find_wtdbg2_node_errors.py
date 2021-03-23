#!/usr/bin/env python3

import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

if len(sys.argv) != 3:
	sys.exit("Need exactly the input file and output prefix.")

nodes_path = sys.argv[1]
output_prefix = sys.argv[2]

print("nodes_path = {}".format(nodes_path))
print("output_prefix = {}".format(output_prefix))

class Align:
	def __init__(self, node_id, read_id, contig_id, read_start, read_end, align_start, align_end, forward):
		self.node_id = node_id
		self.read_id = read_id
		self.contig_id = contig_id
		self.read_start = read_start
		self.read_end = read_end
		self.align_start = align_start
		self.align_end = align_end
		self.forward = forward

with open(nodes_path, 'r') as nodes_file:
	node_align_map = {}
	for line in nodes_file:
		columns = [s for s in line.split("\t") if len(s) > 0]

		node_id = columns[0]

		aligns = []
		for column in columns[2:]:
			column = column.strip()
			#print(column)
			fields = column.split("_")
			coordinate_fields = fields[-5:]
			
			read_id = fields[0]
			contig_id = "_".join(fields[1:-5])

#NC_002695.2_3080719_3105661_F_41_4 -> 3091215_3092239
#NC_002695.2_3110074_3085294_R_69_4 -> 3091386_3092410 (shift 171)
#NC_002695.2_3115177_3090418_R_89_4 -> 3091369_3092393 (shift 154)
#NC_002695.2_3086190_3110151_F_20_4 -> 3091310_3092334 (shift 95)

			read_start = int(coordinate_fields[0])
			read_end = int(coordinate_fields[1])
			forward = coordinate_fields[2] == "F"
			align_read_offset = int(coordinate_fields[3]) * 256

			star = False
			if coordinate_fields[4].endswith("*"):
				coordinate_fields[4] = coordinate_fields[4][:-1]
				star = True
			align_len = int(coordinate_fields[4]) * 256

			if read_start < read_end:
				align_start = read_start + align_read_offset
				align_end = align_start + align_len
			else:
				read_start, read_end = read_end, read_start
				forward = not forward
				align_end = read_end - align_read_offset
				align_start = align_end - align_len

			align = Align(node_id, read_id, contig_id, read_start, read_end, align_start, align_end, forward)
			aligns.append(align)

		node_align_map[node_id] = aligns

correct_nodes = 0
transitively_correct_nodes = 0

deviations = []
for node_id, aligns in node_align_map.items():
	aligns.sort(key = lambda a: (a.align_start, a.align_end))

	min_start = aligns[0].align_start
	max_start = aligns[0].align_start
	min_end = aligns[0].align_end
	max_end = aligns[0].align_end
	contigs = set()
	contigs.add(aligns[0].contig_id)
	is_transitively_correct = True

	for last_align, align in zip(aligns[:-1], aligns[1:]):
		min_start = min(min_start, align.align_start)
		max_start = max(max_start, align.align_start)
		min_end = min(min_end, align.align_end)
		max_end = max(max_end, align.align_end)
		contigs.add(align.contig_id)

		if abs(align.align_start - last_align.align_start) > 128 or abs(align.align_end - last_align.align_end) > 128:
			is_transitively_correct = False

	start_deviation = max_start - min_start
	end_deviation = max_end - min_end
	deviation = max(start_deviation, end_deviation)
	if len(contigs) > 1:
		#print(contigs)
		deviation = None

	if deviation is not None and deviation < 256:
		correct_nodes += 1
	deviations.append(deviation)

	if is_transitively_correct:
		transitively_correct_nodes += 1

print("{} nodes".format(len(node_align_map)))
print("{:.0f}% of nodes are correct".format(correct_nodes * 100.0 / len(node_align_map)))
print("{:.0f}% of nodes are transitively correct".format(transitively_correct_nodes * 100.0 / len(node_align_map)))
print("{:.0f}% of nodes have large errors".format(len([d for d in deviations if d is not None and d >= 1024]) * 100.0 / len(node_align_map)))
print("{:.0f}% of nodes join contigs".format(len([d for d in deviations if d is None]) * 100.0 / len(node_align_map)))

deviations = [d if d > 0 else 1 for d in deviations if d is not None]
deviation_df = pd.DataFrame({"deviation": deviations})

fig, axes = plt.subplots(1, 2)
deviation_histogram = sns.histplot(ax = axes[0], data = deviation_df, x = "deviation", stat = "count", log_scale = 2)
deviation_histogram = sns.histplot(ax = axes[1], data = deviation_df[(deviation_df["deviation"] >= 128) & (deviation_df["deviation"] <= 512)], x = "deviation", stat = "count", log_scale = 2)
fig.savefig(os.path.join(output_prefix, "deviation_histogram.pdf"))