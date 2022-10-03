#!/usr/bin/env python3

import os, sys, json, pathlib, shutil

base_dir = sys.argv[1]

# collect subdirs

existing_key_value_pairs = {}
datas = {}
for subdir, dirs, files in os.walk(base_dir):
    if "report.tex" not in files:
        continue

    data_dir = os.path.relpath(subdir, base_dir)
    if data_dir.startswith("manual") or data_dir.startswith("/manual"):
        continue

    data_str = data_dir.replace("/_/", "")

    try:
        data = json.loads(data_str)
    except json.decoder.JSONDecodeError:
        print(data_str)
        raise

    datas[subdir] = data

    for key, value in data.items():
        existing_key_value_pairs.setdefault(key, set()).add(value)

keys_with_multiple_values = dict([(key, value) for key, value in existing_key_value_pairs.items() if len(value) >= 2])
keys_with_single_value = dict([(key, list(value)[0]) for key, value in existing_key_value_pairs.items() if len(value) == 1])
print(keys_with_multiple_values)

def parse_report(report_file, data):
    report = ''.join(report_file.readlines())
    report_data = {}

    while True:
        index = report.find("\\begin{tabular}")
        if index == -1:
            break
        report = report[index:].strip()

        PARAMETER = "Parameter"
        index = report.find(PARAMETER)
        assert index != -1
        report = report[index:].strip()

        index = report.find("\n")
        assert index != -1
        assembler_line = report[:index]
        report = report[index:].strip()

        assemblers = [assembler.strip() for assembler in assembler_line.split("&")][1:]
        if len(assemblers) == 0:
            continue
        assemblers[-1] = assemblers[-1].replace("\\hline", "").strip().rstrip("\\ ")
        print(assemblers)

        while True:
            index = report.find("\n")
            assert index != -1
            data_line = report[:index]
            report = report[index:].strip()

            if data_line == "\\end{tabular}":
                break
            assert "\\end{tabular}" not in data_line

            if "&" not in data_line:
                continue
            entries = [entry.strip() for entry in data_line.split("&")]
            assert len(entries) == len(assemblers) + 1, f"entries: {entries}\nassemblers: {assemblers}"
            entries[-1] = entries[-1].replace("\\hline", "").strip().rstrip("\\ ")

            line_key = entries[0]
            entries = entries[1:]

            for assembler, entry in zip(assemblers, entries):
                report_data.setdefault(assembler, dict())[line_key] = entry

    return report_data

def generate_report(subdir, data):
    with open(os.path.join(subdir, "report.tex"), 'r') as report_file:
        return parse_report(report_file, data)

def encode_dir(plain):
    STEP = 250
    encoded = []
    for index in range(0, len(plain), STEP):
        if index + STEP <= len(plain):
            encoded.append(plain[index:index+STEP])
        else:
            encoded.append(plain[index:])
    return "/_/".join(encoded)

PARAMETERS = {
    "\\# contigs": "\\# contigs",
    "\\# unique misassemblies": "\\# unique misassemblies",
    "Genome fraction (\\%)": "Genome fraction (\\%)",
    "Duplication ratio": "Duplication ratio",
    "EA50max": "EA50max",
    "EA75max": "EA75max",
    "time [s]": "time [s]",
    "mem [GiB]": "mem [GiB]",
    "hoco time": "hoco time [s]",
    "hodeco time": "hodeco time [s]",
    "hodeco mem": "hodeco mem [GiB]",
    "trivial\\_omnitigs time": "YV time [s]",
    "trivial\\_omnitigs mem": "YV mem [GiB]",
}

#w2 noho & w2 & w2 sfa & w2 YV sfa & w2 YV & w2 frg YV & w2 frg YV fcs=2 & flye & flye sac & flye YV & flye YV sac & hifiasm & mdbg & lja & HiCanu
ASSEMBLERS = {
    "w2": "wtdbg2",
    "w2": "wtdbg2 hoco",
    "w2 frg YV": "wtdbg2 hoco YV",
    "flye": "flye",
    "flye sac": "flye sac",
    "flye YV sac": "flye sac YV",
    "hifiasm": "hifiasm",
    "mdbg": "mdbg",
    "lja": "lja",
    "HiCanu": "hicanu",
}

manual_dir = os.path.join(base_dir, "manual")
shutil.rmtree(manual_dir, ignore_errors=True)

for subdir, data in datas.items():
    report_data = generate_report(subdir, data)

    compressed_variant_dir_dict = {}
    for key in keys_with_multiple_values:
        compressed_variant_dir_dict[key] = data[key]
    compressed_variant_dir = json.dumps(compressed_variant_dir_dict, separators = (',', ':'))
    compressed_variant_dir = encode_dir(compressed_variant_dir)

    manual_subdir = os.path.join(manual_dir, compressed_variant_dir)
    manual_file = os.path.join(manual_subdir, "report.tex")
    print(f"Generating {manual_file}")
    pathlib.Path(manual_subdir).mkdir(parents=True, exist_ok=True)

    with open(manual_file, 'w') as output:
        caption = ", ".join([f"{key}={value}" for key, value in compressed_variant_dir_dict.items()])

        output.write("\\begin{table}[htb]\n")
        output.write("\\centering\n")
        output.write(f"\\caption{{{caption}}}\n")

        output.write("\\begin{tabular}{l")
        for i in range(len(ASSEMBLERS)):
            output.write("r")
        output.write("}\n")

        output.write(" & ")
        output.write(" & ".join(ASSEMBLERS.values()))
        output.write("\\\\\\hline\n")

        for parameter_key, parameter_name in PARAMETERS.items():
            output.write(parameter_name)
            for assembler_key, assembler_name in ASSEMBLERS.items():
                output.write(" & ")
                output.write(report_data[assembler_key][parameter_key])
            output.write("\\\\\n")

        output.write("\\end{tabular}\n")
        output.write("\\end{table}\n")

combined_output_file = os.path.join(manual_dir, "all.tex")
print(f"Generating {combined_output_file}")

with open(combined_output_file, 'w') as output:
    output.write("\\documentclass[10pt,a4paper]{article}\n\\usepackage[margin=0pt]{geometry}\n\\usepackage{lmodern}\n\\usepackage[T1]{fontenc}\n\\usepackage{graphicx}\n")
    output.write("\\begin{document}\n\n")

    output.write("\\textbf{Global parameters}\n")
    output.write("\\begin{itemize}\n")
    for key, value in keys_with_single_value.items():
        if isinstance(key, str):
            key = key.replace("_", "\\_")
        if isinstance(value, str):
            value = value.replace("_", "\\_")
        output.write(f"\\item {key}={value}\n")
    output.write("\\end{itemize}\n")

    for subdir, dirs, files in os.walk(manual_dir):
        if "report.tex" not in files:
            continue

        output.write("\n")
        with open(os.path.join(subdir, "report.tex"), 'r') as report_file:
            for line in report_file.readlines():
                output.write(line.replace("_", "\\_"))

    output.write("\n\\end{document}\n")
