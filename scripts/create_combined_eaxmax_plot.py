#!/usr/bin/python3

import sys

input_shortnames = sys.argv[1:-1:2]
input_quasts = sys.argv[2:-1:2]
output_file = sys.argv[-1]

from os import path


import pandas
df = pandas.DataFrame(columns = ["experiment", "x", "y"])

for shortname, quast_path in zip(input_shortnames, input_quasts):
    file = path.join(quast_path, "aligned_stats/EAxmax_plot.csv")
    frame = pandas.read_csv(file, names=["x", "y"])
    frame["experiment"] = shortname
    df = df.append(frame)

import seaborn as sns
import matplotlib.pyplot as plt
plot = sns.lineplot(data=df, x="x", y="y", hue="experiment")
plt.savefig(output_file)