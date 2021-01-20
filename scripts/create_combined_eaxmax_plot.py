#!/usr/bin/python3

import sys

input_quasts = sys.argv[1:-1]
output_file = sys.argv[-1]

from os import path


import pandas
df = pandas.DataFrame(columns = ["experiment", "x", "y"])

for quast_path in input_quasts:
    name = dir[7:-6]
    file = path.join(quast_path, "aligned_stats/EAxmax_plot.csv")
    frame = pandas.read_csv(file, names=["x", "y"])
    frame["experiment"] = name
    df = df.append(frame)

import seaborn as sns
import matplotlib.pyplot as plt
plot = sns.lineplot(data=df, x="x", y="y", hue="experiment")
plt.savefig(output_file)