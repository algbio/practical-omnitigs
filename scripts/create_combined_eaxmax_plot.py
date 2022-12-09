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

df["EAxmax [million bp]"] = df["EAxmax"] / 1_000_000
try:
    max_eaxmax = max(df["EAxmax"])
except ValueError:
    max_eaxmax = 0
print(df.to_string())

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

legend_visible = True # "D.melanogaster_plot" not in output_file

fig = plt.figure(figsize = (5, 3.5))

SPECIAL_INPUT_SHORTNAMES = set(["F", "F-so", "W", "W-so"])
if set(input_shortnames) == SPECIAL_INPUT_SHORTNAMES:
    palette = {
        "W": "black",
        "W-so": "black",
        "F": "orange",
        "F-so": "orange",
    }
    dashes = {
        "W": "",
        "W-so": (1, 1),
        "F": "",
        "F-so": (1, 1),
    }
    ax = sns.lineplot(data=df, x="x", y="EAxmax [million bp]", hue="Assembler", palette=palette, style="Assembler", dashes=dashes)
else:
    ax = sns.lineplot(data=df, x="x", y="EAxmax [million bp]", hue="Assembler")

if max_eaxmax < 5_000_000:
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: '{:,.1f}'.format(y)))
else:
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: '{:,.0f}'.format(y)))

if ax.get_legend() is not None:
    ax.get_legend().set_visible(legend_visible)

fig.add_axes(ax)
plt.savefig(output_file, bbox_inches="tight")