#!/usr/bin/python3

import sys

input_shortnames = sys.argv[1:-1:2]
input_quast_csvs = sys.argv[2:-1:2]
output_file = sys.argv[-1]

from os import path


import pandas
df = pandas.DataFrame(columns = ["Assembler", "x", "EAxmax"])

for shortname, quast_csv in zip(input_shortnames, input_quast_csvs):
    frame = pandas.read_csv(quast_csv, names=["x", "EAxmax"])
    frame["Assembler"] = shortname
    df = df.append(frame)

df["EAxmax [million bp]"] = df["EAxmax"] / 1000000

import seaborn as sns
import matplotlib.pyplot as plt
plot = sns.lineplot(data=df, x="x", y="EAxmax [million bp]", hue="Assembler")
plt.savefig(output_file)