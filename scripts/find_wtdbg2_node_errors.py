#!/usr/bin/env python3

import sys, os, time
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from matplotlib.lines import Line2D # Imported for legends
import seaborn as sns
import pandas as pd
import scipy.cluster.hierarchy as hierarchy
import numpy as np
import scipy.spatial.distance as ssd
import Levenshtein
from matplotlib.backends.backend_pdf import PdfPages

if len(sys.argv) != 4:
    sys.exit("Need exactly the input nodes file, the input reference file and the output prefix.")

nodes_path = sys.argv[1]
reference_path = sys.argv[2]
output_prefix = sys.argv[3]

print("nodes_path = {}".format(nodes_path))
print("reference_path = {}".format(reference_path))
print("output_prefix = {}".format(output_prefix))

print("Loading reference")
reference_loading_start = time.time()
reference = {}
with open(reference_path, 'r') as reference_file:
    reference_id = None
    sequences = []
    for line in reference_file:
        if line.startswith('>'):
            if reference_id is not None:
                assert len(sequences) > 0, "Found empty reference"
                reference[reference_id] = ''.join(sequences)
                sequences = []

            reference_id = line[1:].strip().split(' ')[0].strip()
            print(reference_id)
        else:
            assert reference_id is not None, "Found sequence without reference id"
            sequences.append(line.strip())

    if reference_id is not None:
        assert len(sequences) > 0, "Found empty reference"
        reference[reference_id] = ''.join(sequences)
    reference_id = None
    sequences = None

print("Reference loading took {}s".format(time.time() - reference_loading_start))

class Align:
    def __init__(self, node_id, read_id, contig_id, read_start, read_end, align_start, align_end, forward, star):
        self.node_id = node_id
        self.read_id = read_id
        self.contig_id = contig_id
        self.read_start = read_start
        self.read_end = read_end
        self.read_len = read_end - read_start
        self.align_start = align_start
        self.align_end = align_end
        self.align_center = (align_start + align_end) // 2
        self.forward = forward
        self.star = star

print("Loading aligns")
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

            align = Align(node_id, read_id, contig_id, read_start, read_end, align_start, align_end, forward, star)
            aligns.append(align)

        node_align_map[node_id] = aligns

print("\n === Investigating nodes ===\n")
correct_nodes = 0
transitively_correct_nodes = 0

deviations = []
large_error_nodes = {}
transitive_subnode_intervals = []
for node_id, aligns in node_align_map.items():
    aligns.sort(key = lambda a: (a.align_center, a.align_start, a.align_end))

    min_start = aligns[0].align_start
    max_start = aligns[0].align_start
    min_end = aligns[0].align_end
    max_end = aligns[0].align_end
    last_transitive_start = min_start
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
            transitive_subnode_intervals.append((last_transitive_start, last_align.align_end))
            last_transitive_start = align.align_start

    transitive_subnode_intervals.append((last_transitive_start, max_end))

    start_deviation = max_start - min_start
    end_deviation = max_end - min_end
    deviation = max(start_deviation, end_deviation)
    if len(contigs) > 1:
        #print(contigs)
        deviation = None

    if deviation is not None and deviation < 256:
        correct_nodes += 1
    deviations.append(deviation)

    if deviation >= 1024:
        large_error_nodes[node_id] = aligns

    if is_transitively_correct:
        transitively_correct_nodes += 1

transitive_subnode_intervals.sort()

print("{} nodes".format(len(node_align_map)))
print("{:.0f}% of nodes are correct".format(correct_nodes * 100.0 / len(node_align_map)))
print("{:.0f}% of nodes are transitively correct".format(transitively_correct_nodes * 100.0 / len(node_align_map)))
print("{:.0f}% of nodes have large errors".format(len([d for d in deviations if d is not None and d >= 1024]) * 100.0 / len(node_align_map)))
print("{:.0f}% of nodes join contigs".format(len([d for d in deviations if d is None]) * 100.0 / len(node_align_map)))

deviations = [d if d > 0 else 1 for d in deviations if d is not None]
deviation_df = pd.DataFrame({"deviation": deviations})

fig, axes = plt.subplots(1, 2, sharey = True)
deviation_histogram = sns.histplot(ax = axes[0], data = deviation_df, x = "deviation", stat = "count", log_scale = 2)
deviation_histogram = sns.histplot(ax = axes[1], data = deviation_df[(deviation_df["deviation"] >= 128) & (deviation_df["deviation"] <= 512)], x = "deviation", stat = "count", log_scale = 2)

path = os.path.join(output_prefix, "deviation_histogram.pdf")
print("Saving {}".format(path))
fig.savefig(path)

print("\n === Investigating large error nodes ===\n")

large_error_dendrograms_path = os.path.join(output_prefix, "large_error_dendrograms.pdf")
large_error_dendrograms_pp = PdfPages(large_error_dendrograms_path)

large_error_figs = [None] * len(large_error_nodes)
large_error_axes = [None] * len(large_error_nodes)
for i in range(len(large_error_nodes)):
    large_error_figs[i], large_error_axes[i] = plt.subplots(2, 1, figsize = (5, 10), constrained_layout = True)

for i, (node_id, aligns) in enumerate(large_error_nodes.items()):
    align_centers_np = np.array([[align.align_center, 0] for align in aligns])
    #print("len(align_centers): {}".format(len(align_centers)))
    linked = hierarchy.linkage(align_centers_np, 'single')

    labelList = range(1, len(align_centers_np) + 1)

    #print("linked.shape: {}".format(linked.shape))
    hierarchy.dendrogram(linked, ax = large_error_axes[i][0],
                orientation='top',
                labels=labelList,
                distance_sort='descending',
                show_leaf_counts=True)
    large_error_axes[i][0].set_ylim([10, 2e6])
    large_error_axes[i][0].set_yscale("log")

# Find transitive subnodes
largest_transitive_subnode_fraction = []
for i, (node_id, aligns) in enumerate(large_error_nodes.items()):
    transitive_subnodes = [[aligns[0]]]
    for align in aligns[1:]:
        if align.align_center - transitive_subnodes[-1][-1].align_center <= 128:
            transitive_subnodes[-1].append(align)
        else:
            transitive_subnodes.append([align])

    transitive_subnode_lens = [len(transitive_subnode) for transitive_subnode in transitive_subnodes]
    transitive_subnode_lens.sort(reverse = True)
    largest_transitive_subnode_fraction.append(transitive_subnode_lens[0] / len(aligns))

largest_transitive_subnode_fraction_df = pd.DataFrame({"largest_transitive_subnode_fraction": largest_transitive_subnode_fraction})
fig, axes = plt.subplots(1, 1, figsize=(5, 5))
largest_transitive_subnode_fraction_histogram = sns.histplot(ax = axes, data = largest_transitive_subnode_fraction_df, x = "largest_transitive_subnode_fraction", stat = "count")
path = os.path.join(output_prefix, "largest_transitive_subnode_fraction_histogram.pdf")
print("Saving {}".format(path))
fig.savefig(path)

# Compute accurate alignments within large error nodes
trans = str.maketrans('ATGC', 'TACG')
def reverse_complement(s):
    return s.translate(trans)

# Get the n largest clusters from the linkage matrix
def get_n_largest_clusters(linkage, n):
    original_id_limit = linkage.shape[0] + 1
    result = [set([linkage.shape[0] + original_id_limit - 1])]
    #print(result)
    #print(linkage)

    for i in reversed(range(linkage.shape[0])):
        row = linkage[i]
        cluster_index = i + original_id_limit
        indices_found = 0
        for result_cluster_index in range(len(result)):
            result_cluster = result[result_cluster_index]
            if cluster_index in result_cluster:
                result_cluster.remove(cluster_index)
                result_cluster.add(int(row[0]))
                if len(result) < n:
                    result.append(set([int(row[1])]))
                else:
                    result_cluster.add(int(row[1]))
                indices_found += 1

        assert indices_found == 1, "Found {} copies of index {}".format(indices_found, cluster_index)
    return result


print("Computing true edit distances")

large_error_nodes_with_transitively_correct_biclustering_count = 0
large_error_nodes_cluster_split_to_transitive_count = []
#fig, axes = plt.subplots(1, len(large_error_nodes), figsize=(5 * len(large_error_nodes), 5), constrained_layout = True, sharey = True)

for plot_index, (node_id, aligns) in enumerate(large_error_nodes.items()):
    print("{:.0f}%".format(plot_index * 100.0 / len(large_error_nodes)))
    distance_matrix = np.zeros((len(aligns), len(aligns)))
    for i, align1 in enumerate(aligns):
        for j, align2 in enumerate(aligns):
            if j >= i:
                continue

            s1 = reference[align1.contig_id][align1.align_start:align1.align_end]
            s2 = reference[align2.contig_id][align2.align_start:align2.align_end]

            if not align1.forward:
                s1 = reverse_complement(s1)
            if not align2.forward:
                s2 = reverse_complement(s2)

            #edit_distance = levenshteinDistance(s1, s2)
            edit_distance = Levenshtein.distance(s1, s2)
            distance_matrix[i][j] = edit_distance
            distance_matrix[j][i] = edit_distance

    #print(distance_matrix)

    # convert the redundant n*n square matrix form into a condensed nC2 array
    distance_matrix = ssd.squareform(distance_matrix) # distance_matrix[{n choose 2}-{n-i choose 2} + (j-i-1)] is the distance between points i and j
    linked = hierarchy.linkage(distance_matrix, 'single')

    labelList = range(1, len(aligns) + 1)

    #print(linked)

    #print("linked.shape: {}".format(linked.shape))
    hierarchy.dendrogram(linked, ax = large_error_axes[plot_index][1],
                orientation='top',
                labels=labelList,
                distance_sort='descending',
                show_leaf_counts=True)
    large_error_axes[plot_index][1].set_ylim([0, 1000])

    is_transitively_correct = False
    split_count = 1
    min_distance = None
    while not is_transitively_correct:
        min_distance = linked[linked.shape[0] - split_count][2]
        split_count += 1
        largest_clusters = get_n_largest_clusters(linked, split_count)
        largest_clusters_aligns = []
        for cluster in largest_clusters:
            cluster_aligns = []
            for index in cluster:
                cluster_aligns.append(aligns[index])
            largest_clusters_aligns.append(cluster_aligns)

        is_transitively_correct = True
        for cluster_aligns in largest_clusters_aligns:
            cluster_aligns.sort(key = lambda align: align.align_center)
            last_center = cluster_aligns[0].align_center
            for align in cluster_aligns[1:]:
                if align.align_center - last_center > 128:
                    is_transitively_correct = False
                    break
                else:
                    last_center = align.align_center


    large_error_nodes_cluster_split_to_transitive_count.append(split_count)

print("Drawing figures for {}".format(large_error_dendrograms_path))
for i in range(len(large_error_nodes)):
    large_error_figs[i].savefig(large_error_dendrograms_pp, format = "pdf")

print("Saving {}".format(large_error_dendrograms_path))
large_error_dendrograms_pp.close()

large_error_nodes_with_transitively_correct_biclustering_count = large_error_nodes_cluster_split_to_transitive_count.count(2)
print("{:.0f}% of large error nodes can be split into two transitive correct nodes by exact hierarchical clustering.".format(large_error_nodes_with_transitively_correct_biclustering_count * 100.0 / len(large_error_nodes)))
large_error_nodes_with_transitively_correct_bi_tri_clustering_count = large_error_nodes_with_transitively_correct_biclustering_count + large_error_nodes_cluster_split_to_transitive_count.count(3)
print("{:.0f}% of large error nodes can be split into two or three transitive correct nodes by exact hierarchical clustering.".format(large_error_nodes_with_transitively_correct_bi_tri_clustering_count * 100.0 / len(large_error_nodes)))
large_error_nodes_with_transitively_correct_bi_tri_quad_clustering_count = large_error_nodes_with_transitively_correct_bi_tri_clustering_count + large_error_nodes_cluster_split_to_transitive_count.count(4)
print("{:.0f}% of large error nodes can be split into two, three or four transitive correct nodes by exact hierarchical clustering.".format(large_error_nodes_with_transitively_correct_bi_tri_quad_clustering_count * 100.0 / len(large_error_nodes)))

large_error_nodes_cluster_split_to_transitive_count_df = pd.DataFrame({"large_error_nodes_cluster_split_to_transitive_count": large_error_nodes_cluster_split_to_transitive_count})
fig, axes = plt.subplots(1, 1, figsize=(5, 5))
large_error_nodes_cluster_split_to_transitive_count_histogram = sns.histplot(ax = axes, data = large_error_nodes_cluster_split_to_transitive_count_df, x = "large_error_nodes_cluster_split_to_transitive_count", stat = "count", binwidth = 1)
path = os.path.join(output_prefix, "large_error_nodes_cluster_split_to_transitive_count_histogram.pdf")
print("Saving {}".format(path))
fig.savefig(path)

# Investigate nodes that should have been glued, but have not been

unglued_nodes = {}
bucket_node_map = {}

def div_ceil(x, y):
    return 1 + ((x - 1) // y)

for contig_id, sequence in reference.items():
    bucket_node_map[contig_id] = []
    for i in range(div_ceil(len(sequence), 128) + 1):
        bucket_node_map[contig_id].append(set())

for node_id, aligns in node_align_map.items():
    for align in aligns:
        for i in range((align.align_start // 128) + 1, div_ceil(align.align_end, 128)):
            bucket_node_map[align.contig_id][i].add(node_id)

bucket_count = 0
empty_bucket_count = 0
underglued_bucket_count = 0
bucket_lens = []

for contig_id, buckets in bucket_node_map.items():
    for bucket_index, bucket in enumerate(buckets):
        bucket_lens.append(len(bucket))
        bucket_count += 1
        if len(bucket) > 1:
            underglued_bucket_count += 1
            #print(bucket)
        elif len(bucket) == 0:
            empty_bucket_count += 1

        if len(bucket) > 10:
            print("contig {} bucket {} has len {}".format(contig_id, bucket_index, len(bucket)))

print("{:.0f}% of reference positions are not part of any node".format(empty_bucket_count * 100.0 / bucket_count))
print("{:.0f}% of reference positions are underglued".format(underglued_bucket_count * 100.0 / bucket_count))

bucket_lens_df = pd.DataFrame({"bucket_lens": bucket_lens})
fig, axes = plt.subplots(1, 1, figsize=(5, 5))
underglued_bucket_histogram = sns.histplot(ax = axes, data = bucket_lens_df, x = "bucket_lens", stat = "count", binwidth = 1)
path = os.path.join(output_prefix, "underglued_bucket_histogram.pdf")
print("Saving {}".format(path))
fig.savefig(path)

# Investigate node align overlaps

intervals = transitive_subnode_intervals[:10]

#viridis = plt.cm.get_cmap('viridis', num_intervals)
#colors = np.array([viridis(idx / num_intervals) for idx in range(len(intervals))])

# Prepare the input data in correct format for LineCollection 
lines = []
layer_ends = [0]
for interval in intervals:
    layer_index = None
    for layer_end_index in range(len(layer_ends)):
        if layer_ends[layer_end_index] <= interval[0] - 1.0:
            layer_index = layer_end_index
            break

    if layer_index is None:
        layer_index = len(layer_ends)
        layer_ends.append(None)

    lines.append([(interval[0], layer_index), (interval[1], layer_index)])
    layer_ends[layer_index] = interval[1]

lc = mc.LineCollection(lines, linewidths = 1)
fig, ax = plt.subplots()
ax.add_collection(lc)
ax.margins(0.1)
#plt.yticks([], [])

# Adding annotations
#for i, x in enumerate(intervals):
#    plt.text(x[0], i+0.1, x[0])
#    plt.text(x[1], i+0.1, x[1])

path = os.path.join(output_prefix, "align_overlaps.pdf")
print("Saving {}".format(path))
plt.savefig(path)