#!/usr/bin/python3

import sys

if len(sys.argv) != 3:
    sys.exit("Wrong number of arguments, expected 2")

input_path = sys.argv[1]
output_file = sys.argv[2]

from os import listdir, path

input_directories = [f for f in listdir(input_path) if path.isdir(path.join(input_path, f)) and f.endswith(".quast")]

import pandas
df = pandas.DataFrame(columns = ["experiment", "x", "y"])

for dir in input_directories:
    name = dir[7:-6]
    path = path.join(input_path, dir, "aligned_stats/EAxmax_plot.csv")
    frame = pandas.read_csv(path, names=["x", "y"])
    frame["experiment"] = name
    df = df.append(frame)

import seaborn as sns
plot = sns.lineplot(data=df, x="x", y="y", hue="experiment")
plot.savefig(output_file)