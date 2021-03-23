import pathlib, itertools, sys

###############################
###### Preprocess Config ######
###############################

print("Preprocessing config", flush = True)

import itertools
import sys
import re
import json
import traceback

configfile: "config/default.yml"

# Allow to configure to use conda from the config file
if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

workflow.global_resources["contigvalidator"] = 1
workflow.global_resources["concorde"] = 1

DATADIR = "data/"
if "datadir" in config:
    DATADIR = config["datadir"]
if DATADIR[-1] != "/":
    DATADIR += "/"

PROGRAMDIR = "data/"
if "programdir" in config:
    PROGRAMDIR = config["programdir"]
if PROGRAMDIR[-1] != "/":
    PROGRAMDIR += "/"

REPORTDIR = "data/reports/"
if "reportdir" in config:
    REPORTDIR = config["reportdir"]
if REPORTDIR[-1] != "/":
    REPORTDIR += "/"

MAX_THREADS = 56
print("Setting MAX_THREADS to " + str(MAX_THREADS), flush = True)

# Preprocess experiments configuration
experiments_bcalm2 = config["experiments"]["bcalm2"]
genomes = config["genomes"]
corrected_genomes = config["corrected_genomes"]
reports = config["reports"]
aggregated_reports = config["aggregated_reports"]

class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'

def safe_format(str, **kwargs):
    return str.format_map(SafeDict(kwargs))

class Arguments(dict):
    @staticmethod
    def from_dict(arguments):
        assert type(arguments) is dict or type(arguments) is OrderedDict

        result = Arguments()
        for key, value in arguments.items():
            if type(value) is dict or type(value) is OrderedDict:
                result[key] = Arguments.from_dict(value)
            else:
                result[key] = value
        return result

    def update(self, other):
        #print("Calling update, other: {}".format(other))

        if other is None:
            return
        assert type(other) is Arguments

        SELECT_PREFIX = "select_"
        for key, value in other.items():
            if key.startswith(SELECT_PREFIX):
                key = key[len(SELECT_PREFIX):]
                select = True
            else:
                select = False

            if key not in self or self[key] is None:
                self[key] = value
            elif type(value) is str or type(value) is int or type(value) is float or type(value) is bool:
                # Selection logic: If only one value is given while we have a dict of values,
                # then select the dict entry with the respective key.
                if type(self[key]) is Arguments:
                    unselect = [k for k in self[key].keys() if k != value]
                    for k in unselect:
                        self[key].pop(k)

                    if len(self[key]) == 0:
                        sys.exit("Unselected all values")
                else:
                    assert type(self[key]) is str or type(self[key]) is int or type(self[key]) is float or type(self[key]) is bool, "type(self[key]) is not str, bool, int or float. key: {}, type(self[key]): {}".format(key, type(self[key]))
                    self[key] = value
            elif type(value) is dict or type(value) is OrderedDict or type(value) is Arguments or value is None:
                assert type(self[key]) is dict or type(self[key]) is OrderedDict or type(self[key]) is Arguments, "wrong type(self[key]): {}".format(type(self[key]))
                # If in selection more, unselect all keys not part of the update.

                unselect = []
                if select:
                    unselect = [k for k in self[key].keys() if k not in value.keys()]
                    for k in unselect:
                        self[key].pop(k)

                self[key].update(value)

                if len(self[key]) == 0 and len(unselect) > 0:
                    sys.exit("Unselected all values: unselect: {}, selected_keys: {}".format(unselect, value.keys()))
            else:
                sys.exit("Cannot merge values of types {} and {}".format(type(self[key]), type(value)))

    def copy(self):
        result = Arguments()
        for key, value in self.items():
            if type(value) is str or type(value) is int or type(value) is float or type(value) is bool or value is None:
                result[key] = value
            elif type(value) is Arguments:
                result[key] = value.copy()
            else:
                sys.exit("Cannot copy value of type {}".format(type(value)))
        return result

    @staticmethod
    def argument_to_argument_string(key, value):
        if value is None or value == True:
            return str(key)
        elif value == False:
            return ""
        else:
            return str(key) + " " + str(value)

    def to_argument_string(self):
        result = ""
        once = True
        for key, value in self.items():
            argument = Arguments.argument_to_argument_string(key, value)
            if argument != "":
                if once:
                    once = False
                else:
                    result += " "
                result += argument
        return result

    def subarguments_name(self, name):
        if name not in self:
            return None

        result = self[name]
        if type(result) is Arguments:
            result = list(result.keys())
            if len(result) != 1:
                return None
            result = result[0]
        return result

    def subarguments_arguments(self, name):
        if name not in self:
            return None

        result = self[name]
        if type(result) is not Arguments:
            return Arguments()

        result = list(result.values())
        if len(result) != 1 or type(result[0]) is not Arguments:
            return Arguments()

        return result[0]

    def read_simulator_name(self):
        return self.subarguments_name("read_simulator")

    def read_simulator_arguments(self):
        return self.subarguments_arguments("read_simulator")

    def assembler_name(self):
        return self.subarguments_name("assembler")

    def assembler_arguments(self):
        return self.subarguments_arguments("assembler")

    def postprocessor_name(self):
        return self.subarguments_name("postprocessor")

    def postprocessor_arguments(self):
        return self.subarguments_arguments("postprocessor")

    def genome(self):
        if "genome" in self:
            return self["genome"]
        else:
            return None

    def __str__(self):
        if self is None:
            return "{}"

        string = json.dumps(self, separators = (',', ':'))
        n = 200
        result = "/_/".join([string[i:i+n] for i in range(0, len(string), n)])
        return result

    def from_str(string):
        if string is None:
            return None
        elif string == "None":
            return None

        string = string.replace("/_/", "")
        return Arguments.from_dict(json.loads(string))

    def shortstr(self):
        string = ""
        shortened = self.copy()
        shortened.pop("genome", None)
        cli_arguments = shortened.get("cli_arguments", None)
        if type(cli_arguments) is Arguments:
            for k, v in cli_arguments.items():
                shortened[k] = v
            shortened.pop("cli_arguments")

        for key, value in shortened.items():
            if len(string) > 0:
                string += ","
            string += key.replace("-", "")[0] + ":"

            if type(value) is Arguments:
                string += "(" + value.shortstr() + ")"
            elif type(value) is int or type(value) is float:
                string += str(value)
            elif type(value) is bool:
                assert value == True
                string = string[:-1]
            else:
                string += str(value)[0:2]
        return string

    def retain_raw_assembly_arguments(self):
        if "postprocessor" in self:
            self.pop("postprocessor")

class Algorithm:
    def __init__(self, assembler, arguments):
        self.assembler = assembler
        self.arguments = arguments

class Column:
    def __init__(self, root_arguments, additional_arguments):
        #print("Creating column\nroot: {}\nadditional: {}".format(root_arguments, additional_arguments))
        assert type(root_arguments) is Arguments
        assert type(additional_arguments) is Arguments
        self.arguments = root_arguments.copy()
        self.arguments.update(additional_arguments.copy())
        #print("Result arguments: {}".format(self.arguments))
        self.shortname = self.arguments.pop("shortname")
        self.assembler = self.arguments.assembler_name()
        self.assembler_arguments = self.arguments.assembler_arguments()
        self.genome = self.arguments["genome"]

    def __str__(self):
        return str(self.arguments)

    def from_str(string):
        return Column(Arguments.from_str(string), Arguments())

class ReportFile:
    def __init__(self, arguments, columns):
        self.arguments = arguments
        self.columns = columns
        self.name = str(arguments)
        self.shortname = arguments.shortstr()

    def __str__(self):
        return "ReportFile(name: {}, arguments: {}, columns: {})".format(self.name, self.arguments, self.columns)

class ArgumentMatrix:
    def __init__(self, argument_matrix):
        self.argument_matrix = argument_matrix
        self.len = ArgumentMatrix._compute_length(argument_matrix)

    def __len__(self):
        return self.len

    def _compute_length(argument_matrix):
        if argument_matrix is None:
            return 1

        if type(argument_matrix) is ArgumentMatrix:
            argument_matrix = argument_matrix.argument_matrix

        if type(argument_matrix) is dict or type(argument_matrix) is OrderedDict:
            length = 1
            for key, value in argument_matrix.items():
                length *= ArgumentMatrix._compute_length(value)
            return length

        if type(argument_matrix) is list:
            length = 0
            for value in argument_matrix:
                length += ArgumentMatrix._compute_length(value)
            return length

        if type(argument_matrix) is str or type(argument_matrix) is int or type(argument_matrix) is float or type(argument_matrix) is bool:
            return 1

        sys.exit("Illegal type of argument matrix: {}.\nAllowed are dict, OrderedDict, list, str, bool, int and float".format(type(argument_matrix)))

    def __iter__(self):
        return ArgumentMatrix._iter(self.argument_matrix)

    def _iter(argument_matrix):
        if argument_matrix is None:
            return None

        if type(argument_matrix) is ArgumentMatrix:
            argument_matrix = argument_matrix.argument_matrix

        if type(argument_matrix) is str or type(argument_matrix) is int or type(argument_matrix) is float or type(argument_matrix) is bool:
            #print("Found str, int or float: {}".format(argument_matrix))
            yield argument_matrix

        elif type(argument_matrix) is dict or type(argument_matrix) is OrderedDict or type(argument_matrix) is Arguments:
            #print("Found dict: {}".format(argument_matrix))
            keys = [None] * len(argument_matrix)
            values = [None] * len(argument_matrix)
            for index, (key, value) in enumerate(argument_matrix.items()):
                keys[index] = key
                if value is None:
                    values[index] = [None]
                else:
                    values[index] = list(ArgumentMatrix._iter(value))

            #print("keys: {}".format(keys))
            #print("values: {}".format(values))

            for value_combination in itertools.product(*values):
                result = Arguments()
                for key, value in zip(keys, value_combination):
                    if value is not None:
                        result[key] = value

                if len(result) == 1:
                    if list(result.values())[0] is None:
                        yield list(result.keys())[0]
                        continue

                yield result

        elif type(argument_matrix) is list:
            #print("Found list: {}".format(argument_matrix))
            for value in argument_matrix:
                if value is None:
                    yield None
                else:
                    for list_item in ArgumentMatrix._iter(value):
                        yield list_item

        else:
            sys.exit("Illegal type of argument matrix: {}.\nAllowed are dict, OrderedDict, list, str, bool, int and float".format(type(argument_matrix)))


for report_name, report_definition in reports.items():
    argument_matrix = ArgumentMatrix(report_definition.setdefault("argument_matrix", {}))
    report_definition["argument_matrix"] = argument_matrix

    # print("Matrix of {} has length {}".format(report_name, len(argument_matrix)))
    # entries = list(iter(argument_matrix))
    # print("Entries in matrix:")
    # for entry in entries:
    #     print(entry)

    for arguments in argument_matrix:
        columns = []
        for column_definition in report_definition["columns"]:
            columns.append(Column(arguments, Arguments.from_dict(column_definition)))
        report_file = ReportFile(arguments, columns)
        report_definition.setdefault("report_files", {})[report_file.name] = report_file


import pprint
pp = pprint.PrettyPrinter(indent = 2, width = 200, compact = True)
pp.pprint(reports)
#pp.pprint(aggregated_reports)

#print("Exiting for debugging")
#sys.exit(0)

tests = config["tests"]

for experiment, config in itertools.chain(experiments_bcalm2.items(), tests.items()):
    if not "url" in config:
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        url += experiment.split("_")[0] + "/"
        number = experiment.split("_")[1].split(".")[0]
        url += number[0:3] + "/"
        url += number[3:6] + "/"
        url += number[6:9] + "/"
        url += experiment + "/"
        url += experiment + "_genomic.fna.gz"
        config["url"] = url

# Collect all rust sources
rust_sources = list(map(str, itertools.chain(pathlib.Path('implementation').glob('**/Cargo.toml'), pathlib.Path('implementation').glob('**/*.rs'))))

import datetime
today = datetime.date.today().isoformat()

print("Finished config preprocessing", flush = True)

#########################
###### Directories ######
#########################

GENOME_READS_FORMAT = DATADIR + "genomes/{genome}/reads/reads.fa"
GENOME_SIMULATED_READS_FORMAT = DATADIR + "genomes/{genome}/simulated_reads/{read_simulator_name}/s/{read_simulator_arguments}/reads.fa"
CORRECTION_SHORT_READS_FORMAT = DATADIR + "corrected_reads/{corrected_genome}/reads/correction_short_reads.fa"
GENOME_REFERENCE_FORMAT = DATADIR + "genomes/{genome}/reference/reference.fa"

ALGORITHM_PREFIX_FORMAT = DATADIR + "algorithms/{arguments}/"
WTDBG2_PREFIX_FORMAT = ALGORITHM_PREFIX_FORMAT + "wtdbg2/"
HIFIASM_PREFIX_FORMAT = os.path.join(ALGORITHM_PREFIX_FORMAT, "hifiasm")

QUAST_PREFIX_FORMAT = ALGORITHM_PREFIX_FORMAT + "quast/"

REPORT_PREFIX_FORMAT = REPORTDIR + "{report_name}/s/{report_file_name}/"
AGGREGATED_REPORT_PREFIX_FORMAT = REPORTDIR + "{aggregated_report_name}/"

#################################
###### Global report rules ######
#################################

def get_all_report_files():
    try:
        result = []
        for report_name, report_definition in reports.items():
            for report_file_name in report_definition["report_files"].keys():
                result.append(REPORT_PREFIX_FORMAT.format(report_name = report_name, report_file_name = report_file_name) + "report.pdf")


        for aggregated_report_name in aggregated_reports.keys():
            result.append(AGGREGATED_REPORT_PREFIX_FORMAT.format(aggregated_report_name = aggregated_report_name) + "aggregated-report.pdf")
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_all_report_files_without_human():
    try:
        result = []
        reports_with_human = set()
        for report_name, report_definition in reports.items():
            for report_file_name in report_definition["report_files"].keys():
                append = True
                for column in report_definition["report_files"][report_file_name].columns:
                    if column.arguments.get("genome", None) in ["HG002", "HG002_HiFi", "HG002_ratatosk"]:
                        append = False
                        reports_with_human.add(report_name)
                        break

                if append:
                    result.append(REPORT_PREFIX_FORMAT.format(report_name = report_name, report_file_name = report_file_name) + "report.pdf")

        for aggregated_report_name in aggregated_reports.keys():
            append = False
            pops = []
            for report_name in aggregated_reports[aggregated_report_name]["reports"]:
                if report_name in reports_with_human:
                    pops.append(report_name)
                else:
                    append = True

            for pop in pops:
                aggregated_reports[aggregated_report_name]["reports"].remove(pop)

            if append:
                result.append(AGGREGATED_REPORT_PREFIX_FORMAT.format(aggregated_report_name = aggregated_report_name) + "aggregated-report.pdf")
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule report_all_without_human:
    input:  get_all_report_files_without_human(),
    threads: 1
    resources: mail_type = "END,FAIL,INVALID_DEPEND,REQUEUE"

rule report_all:
    input:  get_all_report_files(),
    threads: 1
    resources: mail_type = "END,FAIL,INVALID_DEPEND,REQUEUE"

###############################
###### Report Generation ######
###############################

localrules: create_single_bcalm2_report_tex
rule create_single_bcalm2_report_tex:
    input:  genome_name = DATADIR + "{dir}/name.txt",
            unitigs = DATADIR + "{dir}/{file}.unitigs.tex",
            unitigs_contigvalidator = DATADIR + "{dir}/{file}.unitigs.contigvalidator",
            unitigs_quast = DATADIR + "{dir}/{file}.unitigs.quast",
            omnitigs = DATADIR + "{dir}/{file}.omnitigs.tex",
            omnitigs_contigvalidator = DATADIR + "{dir}/{file}.omnitigs.contigvalidator",
            omnitigs_quast = DATADIR + "{dir}/{file}.omnitigs.quast",
            trivialomnitigs = DATADIR + "{dir}/{file}.trivialomnitigs.tex",
            trivialomnitigs_contigvalidator = DATADIR + "{dir}/{file}.trivialomnitigs.contigvalidator",
            trivialomnitigs_quast = DATADIR + "{dir}/{file}.trivialomnitigs.quast",
            graphstatistics = DATADIR + "{dir}/{file}.bcalm2.graphstatistics",
            bcalm2_bandage = DATADIR + "{dir}/{file}.bcalm2.bandage.png",
            script = "scripts/convert_validation_outputs_to_latex.py",
    output: DATADIR + "{dir}/{file}.report.tex",
    params: prefix = DATADIR + "{dir}/{file}",
            datadir = DATADIR,
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    shell: "python3 scripts/convert_validation_outputs_to_latex.py '{params.datadir}' '{input.genome_name}' '{input.graphstatistics}' '{input.bcalm2_bandage}' 'none' '{output}' uni '{params.prefix}.unitigs' 'Y-to-V' '{params.prefix}.trivialomnitigs' omni '{params.prefix}.omnitigs'"

### Create single report ###

def get_report_file(report_name, report_file_name):
    try:
        if report_name not in reports:
            raise Exception("report_name not in reports: {}".format(report_name))
        report_files = reports[report_name]["report_files"]
        if report_file_name not in report_files:
            raise Exception("report_file_name not in report_files: {}".format(report_file_name))

        return report_files[report_file_name]
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_report_file_quasts(report_name, report_file_name):
    try:
        report_file = get_report_file(report_name, report_file_name)
        quasts = []
        for column in report_file.columns:
            quasts.append(QUAST_PREFIX_FORMAT.format(arguments = column.arguments))
        return quasts
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_report_file_quasts_from_wildcards(wildcards):
    try:
        return get_report_file_quasts(wildcards.report_name, wildcards.report_file_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_report_file_column_shortnames(report_name, report_file_name):
    try:
        report_file = get_report_file(report_name, report_file_name)
        shortnames = []
        for column in report_file.columns:
            shortnames.append(column.shortname)
        return shortnames
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_report_file_column_shortnames_from_wildcards(wildcards):
    try:
        return get_report_file_column_shortnames(wildcards.report_name, wildcards.report_file_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_report_genome_names(report_name, report_file_name):
    try:
        report_file = get_report_file(report_name, report_file_name)
        genome_names = set()
        for column in report_file.columns:
            genome = column.genome
            genome_names.add(genome)

        return genome_names
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_report_genome_names_from_wildcards(wildcards):
    try:
        return get_report_genome_names(wildcards.report_name, wildcards.report_file_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_single_report_script_column_arguments(report_name, report_file_name):
    try:
        report_file = get_report_file(report_name, report_file_name)
        result = ""
        once = True
        for column in report_file.columns:
            if once:
                once = False
            else:
                result += " "

            arguments = column.arguments
            result += "'" + column.shortname + "' '" + ALGORITHM_PREFIX_FORMAT.format(arguments = arguments) + "'"
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_single_report_script_column_arguments_from_wildcards(wildcards):
    try:
        return get_single_report_script_column_arguments(wildcards.report_name, wildcards.report_file_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: create_single_report_tex
rule create_single_report_tex:
    input:  quasts = get_report_file_quasts_from_wildcards,
            combined_eaxmax_plot = REPORT_PREFIX_FORMAT + "combined_eaxmax_plot.pdf",
            script = "scripts/convert_validation_outputs_to_latex.py",
    output: report = REPORT_PREFIX_FORMAT + "report.tex",
    params: genome_name = lambda wildcards: ", ".join(get_report_genome_names_from_wildcards(wildcards)),
            script_column_arguments = get_single_report_script_column_arguments_from_wildcards,
            name_file = REPORT_PREFIX_FORMAT + "name.txt",
            hashdir = REPORTDIR + "hashdir",
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    shell: """
        echo '{wildcards.report_name} {params.genome_name} {wildcards.report_file_name}' > '{params.name_file}'
        python3 '{input.script}' '{params.hashdir}' '{params.name_file}' 'none' 'none' '{input.combined_eaxmax_plot}' '{output}' {params.script_column_arguments}
        """

### Create aggregated report ###

def get_aggregated_report_file_maps(aggregated_report_name):
    try:
        result = {}
        for report_name in aggregated_reports[aggregated_report_name]["reports"]:
            result.setdefault(report_name, {}).update(reports[report_name]["report_files"])
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def iterate_aggregated_report_file_source_reports(aggregated_report_name):
    try:
        aggregated_report_files = get_aggregated_report_file_maps(aggregated_report_name)
        for report_name, report_file_map in aggregated_report_files.items():
            for report_file_name, report_file_definition in report_file_map.items():
                yield (report_name, report_file_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def iterate_aggregated_report_file_source_reports_short_names(aggregated_report_name):
    try:
        aggregated_report_files = get_aggregated_report_file_maps(aggregated_report_name)
        for report_name, report_file_map in aggregated_report_files.items():
            for report_file_name, report_file_definition in report_file_map.items():
                yield (report_name, report_file_definition.shortname)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_aggregated_report_file_source_report_paths(aggregated_report_name):
    try:
        result = []
        for report_name, report_file_name in iterate_aggregated_report_file_source_reports(aggregated_report_name):
            result.append(REPORT_PREFIX_FORMAT.format(report_name = report_name, report_file_name = report_file_name) + "report.tex")
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_aggregated_report_file_source_report_paths_from_wildcards(wildcards):
    try:
        return get_aggregated_report_file_source_report_paths(wildcards.aggregated_report_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: create_aggregated_report_tex
rule create_aggregated_report_tex:
    input: source_reports = get_aggregated_report_file_source_report_paths_from_wildcards,
           script = "scripts/create_aggregated_wtdbg2_report.py",
    output: file = AGGREGATED_REPORT_PREFIX_FORMAT + "aggregated-report.tex"
    params: source_reports_arg = lambda wildcards: "' '".join(get_aggregated_report_file_source_report_paths_from_wildcards(wildcards)),
            source_report_names_arg = lambda wildcards: "' '".join([report_name + "/" + report_file_name for report_name, report_file_name in iterate_aggregated_report_file_source_reports_short_names(wildcards.aggregated_report_name)]),
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    shell: """
        python3 '{input.script}' --source-reports '{params.source_reports_arg}' --source-report-names '{params.source_report_names_arg}' --output '{output.file}'
        """

localrules: create_combined_eaxmax_graph
rule create_combined_eaxmax_graph:
    input:  quast_csvs = lambda wildcards: [q + "aligned_stats/EAxmax_plot.csv" for q in get_report_file_quasts_from_wildcards(wildcards)],
            script = "scripts/create_combined_eaxmax_plot.py",
    output: REPORT_PREFIX_FORMAT + "combined_eaxmax_plot.pdf",
    params: input_quast_csvs = lambda wildcards, input: "' '".join([shortname + "' '" + quast for shortname, quast in zip(get_report_file_column_shortnames_from_wildcards(wildcards), input.quast_csvs)])
    conda:  "config/conda-seaborn-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output}')"
        python3 '{input.script}' '{params.input_quast_csvs}' '{output}'
        """

localrules: png_to_pdf
rule png_to_pdf:
    input: "{file}.png"
    output: "{file}.image.pdf"
    conda: "config/conda-imagemagick-env.yml"
    threads: 1
    shell: "convert {input} {output}"

localrules: latex
rule latex:
    input: "{subpath}report.tex"
    output: "{subpath}report.pdf"
    conda: "config/conda-latex-env.yml"
    threads: 1
    shell: """
        tectonic '{input}'
        """

rule find_wtdbg2_node_errors:
    input:  nodes = os.path.join(WTDBG2_PREFIX_FORMAT, "wtdbg2.1.nodes"),
            script = "scripts/find_wtdbg2_node_errors.py",
    log:    log = os.path.join(ALGORITHM_PREFIX_FORMAT, "wtdbg2_node_errors", "wtdbg2_node_errors.log"),
    params: output_prefix = os.path.join(ALGORITHM_PREFIX_FORMAT, "wtdbg2_node_errors") + "/",
    conda:  "config/conda-seaborn-env.yml"
    threads: 1
    shell:  "'{input.script}' '{input.nodes}' '{params.output_prefix}' 2>&1 | tee '{log.log}'"

########################
###### Algorithms ######
########################

### bcalm2 ###

rule compute_omnitigs:
    input: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = PROGRAMDIR + "target/release/cli"
    output: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.omnitigs.fa", log = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.omnitigs.fa.log", latex = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.omnitigs.tex"
    threads: 1
    shell: "'{input.binary}' compute-omnitigs --input '{input.file}' --kmer-size {wildcards.k} --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

rule compute_trivial_omnitigs:
    input: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = PROGRAMDIR + "target/release/cli"
    output: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.trivialomnitigs.fa", log = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.trivialomnitigs.fa.log", latex = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.trivialomnitigs.tex"
    threads: 1
    shell: "'{input.binary}' compute-trivial-omnitigs --input '{input.file}' --kmer-size {wildcards.k} --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

rule compute_unitigs:
    input: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = PROGRAMDIR + "target/release/cli"
    output: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.unitigs.fa", log = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.unitigs.fa.log", latex = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.unitigs.tex"
    threads: 1
    shell: "'{input.binary}' compute-unitigs --input '{input.file}' --kmer-size {wildcards.k} --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

### long reads ###

def get_injectable_contigs_rust_cli_command_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        injections = assembler_arguments.get("injections", {})

        if "wtdbg2-inject-unitigs" in injections:
            return "compute-unitigs"
        elif "wtdbg2-inject-trivial-omnitigs" in injections:
            return "compute-trivial-omnitigs --non-scc"
        elif "wtdbg2-inject-omnitigs" in injections:
            return "compute-omnitigs"
        else:
            raise Exception("Missing injection command in assembler arguments: " + str(wildcards))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_genome_reference_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        return GENOME_REFERENCE_FORMAT.format(genome = arguments["genome"])
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_genome_reads_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)

        if arguments.read_simulator_name() is None or arguments.read_simulator_name() == "none":
            return GENOME_READS_FORMAT.format(genome = arguments["genome"])
        else:
            return GENOME_SIMULATED_READS_FORMAT.format(genome = arguments["genome"], read_simulator_name = arguments.read_simulator_name(), read_simulator_arguments = arguments.read_simulator_arguments())
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule compute_injectable_contigs_wtdbg2:
    input:  nodes = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.3.nodes",
            reads = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.3.reads",
            dot = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.3.dot",
            raw_reads = get_genome_reads_from_wildcards,
            binary = PROGRAMDIR + "target/release/cli",
    output: file = ALGORITHM_PREFIX_FORMAT + "injectable_contigs/contigwalks",
            log = ALGORITHM_PREFIX_FORMAT + "injectable_contigs/compute_injectable_contigs.log",
            latex = ALGORITHM_PREFIX_FORMAT + "injectable_contigs/compute_injectable_contigs.tex",
            completed = touch(ALGORITHM_PREFIX_FORMAT + "injectable_contigs/contigwalks.completed"),
    params: command = get_injectable_contigs_rust_cli_command_from_wildcards
    threads: 1
    resources: mem_mb = 48000
    shell: "'{input.binary}' {params.command} --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

def get_injectable_fragment_contigs_rust_cli_command_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        fragment_injections = assembler_arguments.get("fragment_injections", {})

        if "wtdbg2-inject-fragment-unitigs" in fragment_injections:
            return "compute-unitigs"
        elif "wtdbg2-inject-fragment-trivial-omnitigs" in fragment_injections:
            return "compute-trivial-omnitigs --non-scc"
        elif "wtdbg2-inject-fragment-omnitigs" in fragment_injections:
            return "compute-omnitigs"
        else:
            raise Exception("Missing injection command in assembler arguments: " + str(wildcards))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_injectable_fragment_contigs_input_dot_file_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        fragment_injections = assembler_arguments.get("fragment_injections", {})

        result = fragment_injections.get("wtdbg2-inject-fragment-unitigs", None)
        result = fragment_injections.get("wtdbg2-inject-fragment-trivial-omnitigs", result)
        result = fragment_injections.get("wtdbg2-inject-fragment-omnitigs", result)

        if result is None:
            raise Exception("Missing injection command in assembler arguments: " + str(wildcards))

        if result == True:
            return get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.ctg.dot"
        else:
            return get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2." + str(result) + ".frg.dot"
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule compute_injectable_fragment_contigs_wtdbg2:
    input:  dot = get_injectable_fragment_contigs_input_dot_file_from_wildcards,
            binary = PROGRAMDIR + "target/release/cli",
    output: file = ALGORITHM_PREFIX_FORMAT + "injectable_fragment_contigs/contigwalks",
            latex = ALGORITHM_PREFIX_FORMAT + "injectable_fragment_contigs/compute_injectable_contigs.tex",
            completed = touch(ALGORITHM_PREFIX_FORMAT + "injectable_fragment_contigs/contigwalks.completed"),
    log:    log = ALGORITHM_PREFIX_FORMAT + "injectable_fragment_contigs/compute_injectable_contigs.log",
    params: command = get_injectable_fragment_contigs_rust_cli_command_from_wildcards
    threads: 1
    resources: mem_mb = 48000
    shell: "'{input.binary}' {params.command} --file-format dot --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{log.log}'"

#################################
###### Postprocess Contigs ######
#################################

def get_raw_assembly_file_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_name = arguments.assembler_name()

        if assembler_name == "flye":
            result = ALGORITHM_PREFIX_FORMAT + "flye/assembly.fasta"
        elif assembler_name == "wtdbg2":
            result = WTDBG2_PREFIX_FORMAT + "wtdbg2.raw.fa"
        elif assembler_name == "hifiasm":
            result = os.path.join(HIFIASM_PREFIX_FORMAT, "assembly.p_ctg.fa")
        elif assembler_name == "reference":
            result = GENOME_REFERENCE_FORMAT
        else:
            raise Exception("Unknown assembler {} in arguments: {}".format(assembler_name, arguments))

        arguments.retain_raw_assembly_arguments()
        return result.format(genome = arguments.genome(), arguments = str(arguments))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: select_assembler
rule select_assembler:
    input:  raw_assembly_from_assembler = get_raw_assembly_file_from_wildcards,
    output: raw_assembly = ALGORITHM_PREFIX_FORMAT + "raw_assembly.fa",
    threads: 1
    shell: "ln -sr '{input.raw_assembly_from_assembler}' '{output.raw_assembly}'"

def get_assembly_postprocessing_target_file_from_wildcards(wildcards):
    try:
        try:
            arguments = Arguments.from_str(wildcards.arguments)
        except json.decoder.JSONDecodeError:
            traceback.print_exc()
            raise Exception("JSONDecodeError: {}".format(wildcards.arguments))

        postprocessor_name = arguments.postprocessor_name()
        if postprocessor_name == "contigbreaker":
            return ALGORITHM_PREFIX_FORMAT + "contigbreaker/broken_contigs.fa"
        else:
            return ALGORITHM_PREFIX_FORMAT + "raw_assembly.fa"
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: request_assembly_postprocessing
rule request_assembly_postprocessing:
    input: get_assembly_postprocessing_target_file_from_wildcards
    output: assembly = ALGORITHM_PREFIX_FORMAT + "assembly.fa"
    threads: 1
    shell: "ln -sr '{input}' '{output}'"

def get_source_genome_properties_from_wildcards(wildcards):
    try:
        if hasattr(wildcards, "genome"):
            genome_name = wildcards.genome
        elif hasattr(wildcards, "corrected_genome"):
            genome_name = wildcards.corrected_genome
        elif hasattr(wildcards, "arguments"):
            arguments = Arguments.from_str(wildcards.arguments)
            genome_name = arguments.genome()

            if genome_name is None:
                raise Exception("Arguments has no 'genome' attribute: {}".format(arguments))
        else:
            raise Exception("Wildcards has no 'genome', 'corrected_genome' or 'arguments' attribute")

        if genome_name in genomes:
            pass
        elif corrected_genomes is not None and genome_name in corrected_genomes:
            genome_name = corrected_genomes[genome_name]["source_genome"]
        else:
            raise Exception("Genome name not found: " + str(genome_name))
        return genomes[genome_name]
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def compute_genome_mem_mb_from_wildcards(wildcards, base_mem_mb):
    try:
        genome_properties = get_source_genome_properties_from_wildcards(wildcards)

        if "assembly_mem_factor" in genome_properties:
            return int(float(genome_properties["assembly_mem_factor"]) * float(base_mem_mb))
        else:
            return base_mem_mb
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def compute_genome_time_min_from_wildcards(wildcards, base_time_min):
    try:
        genome_properties = get_source_genome_properties_from_wildcards(wildcards)

        if "assembly_time_factor" in genome_properties:
            return int(float(genome_properties["assembly_time_factor"]) * float(base_time_min))
        else:
            return base_time_min
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def compute_genome_queue_from_wildcards(wildcards, base_time_min, base_mem_mb = 0):
    try:
        time = compute_genome_time_min_from_wildcards(wildcards, base_time_min)
        mem = compute_genome_mem_mb_from_wildcards(wildcards, base_mem_mb)

        if mem >= 250000:
            return "bigmem"
        elif time <= 1440:
            return "short"
        elif time <= 1440 * 14:
            return "long"
        elif time <= 1440 * 60:
            return "extralong"
        else:
            sys.exit("No applicable queue for runtime " + str(time) + " (wildcards: " + str(wildcards) + ")")
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule run_contigbreaker:
    input:  contigs = ALGORITHM_PREFIX_FORMAT + "raw_assembly.fa",
            reads = get_genome_reads_from_wildcards,
            script = "tools/contigbreaker/contigbreaker.py"
    output: broken_contigs = ALGORITHM_PREFIX_FORMAT + "contigbreaker/broken_contigs.fa",
            completed = touch(ALGORITHM_PREFIX_FORMAT + "contigbreaker/broken_contigs.fa.completed"),
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 60000),
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               cpus = MAX_THREADS,
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720),
    conda: "tools/contigbreaker/environment.yml"
    shell: "'{input.script}' --threads {threads} --input-contigs '{input.contigs}' --input-reads '{input.reads}' --output-contigs '{output.broken_contigs}'"

####################
###### wtdbg2 ######
####################

localrules: install_wtdbg2
rule install_wtdbg2:
    output: kbm2 = "external-software/wtdbg2/kbm2", pgzf = "external-software/wtdbg2/pgzf", wtdbg2 = "external-software/wtdbg2/wtdbg2", wtdbg_cns = "external-software/wtdbg2/wtdbg-cns", wtpoa_cns = "external-software/wtdbg2/wtpoa-cns"
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
    mkdir -p external-software
    cd external-software

    git clone https://github.com/sebschmi/wtdbg2.git
    cd wtdbg2
    git checkout c8403f562f3b999bb514ba3e9020007bcf01391c
    make
    """

def get_wtdbg2_args_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        cli_arguments = assembler_arguments.get("cli_arguments", None)
        if cli_arguments is None:
            return ""

        return cli_arguments.to_argument_string()
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_genome_len_from_wildcards(wildcards):
    try:
        genome_properties = get_source_genome_properties_from_wildcards(wildcards)
        return str(genome_properties["genome_length"])
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_caching_arguments_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        cli_arguments = assembler_arguments.get("cli_arguments", {})

        cli_arguments.pop("--skip-fragment-assembly", None)
        cli_arguments.pop("--fragment-correction-steps", None)
        assembler_arguments.pop("injections", None)
        assembler_arguments.pop("fragment_injections", None)
        return arguments
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_caching_prefix_from_wildcards(wildcards):
    try:
        arguments = get_wtdbg2_caching_arguments_from_wildcards(wildcards)
        return WTDBG2_PREFIX_FORMAT.format(arguments = str(arguments))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_injectable_contigs_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        injections = assembler_arguments.get("injections", {})

        if "wtdbg2-inject-unitigs" in injections or "wtdbg2-inject-trivial-omnitigs" in injections or "wtdbg2-inject-omnitigs" in injections:
            return ALGORITHM_PREFIX_FORMAT + "injectable_contigs/contigwalks"

        return []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def is_wtdbg2_injecting_contigs_from_wildcards(wildcards):
    try:
        return get_wtdbg2_injectable_contigs_from_wildcards(wildcards) != []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_injectable_fragment_contigs_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        injections = assembler_arguments.get("fragment_injections", {})

        if "wtdbg2-inject-fragment-unitigs" in injections or "wtdbg2-inject-fragment-trivial-omnitigs" in injections or "wtdbg2-inject-fragment-omnitigs" in injections:
            return ALGORITHM_PREFIX_FORMAT + "injectable_fragment_contigs/contigwalks"

        return []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def is_wtdbg2_injecting_fragment_contigs_from_wildcards(wildcards):
    try:
        return get_wtdbg2_injectable_fragment_contigs_from_wildcards(wildcards) != []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def is_wtdbg2_using_cache_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        cli_arguments = assembler_arguments.get("cli_arguments", {})
        injections = assembler_arguments.get("injections", {})
        fragment_injections = assembler_arguments.get("fragment_injections", {})

        return len(injections) > 0 or len(fragment_injections) > 0 or "--skip-fragment-assembly" in cli_arguments
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_cached_nodes_from_wildcards(wildcards):
    try:
        if is_wtdbg2_using_cache_from_wildcards(wildcards):
            return get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.1.nodes"

        return []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_cached_clips_from_wildcards(wildcards):
    try:
        if is_wtdbg2_using_cache_from_wildcards(wildcards):
            return get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.clps"

        return []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_wtdbg2_cached_kbm_from_wildcards(wildcards):
    try:
        if is_wtdbg2_using_cache_from_wildcards(wildcards):
            return get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.kbm"

        return []
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule wtdbg2:
    input:  reads = get_genome_reads_from_wildcards,
            contigs = get_wtdbg2_injectable_contigs_from_wildcards,
            fragment_contigs = get_wtdbg2_injectable_fragment_contigs_from_wildcards,
            cached_nodes = get_wtdbg2_cached_nodes_from_wildcards,
            cached_clips = get_wtdbg2_cached_clips_from_wildcards,
            cached_kbm = get_wtdbg2_cached_kbm_from_wildcards,
            binary = "external-software/wtdbg2/wtdbg2",
    output: original_nodes = WTDBG2_PREFIX_FORMAT + "wtdbg2.1.nodes",
            nodes = WTDBG2_PREFIX_FORMAT + "wtdbg2.3.nodes",
            reads = WTDBG2_PREFIX_FORMAT + "wtdbg2.3.reads",
            dot = WTDBG2_PREFIX_FORMAT + "wtdbg2.3.dot.gz",
            ctg_dot = WTDBG2_PREFIX_FORMAT + "wtdbg2.ctg.dot.gz",
            clips = WTDBG2_PREFIX_FORMAT + "wtdbg2.clps",
            kbm = WTDBG2_PREFIX_FORMAT + "wtdbg2.kbm",
            ctg_lay = WTDBG2_PREFIX_FORMAT + "wtdbg2.ctg.lay.gz",
            frg_dot = [WTDBG2_PREFIX_FORMAT + "wtdbg2.{}.frg.dot.gz".format(i) for i in range(1, 11)],
            log = WTDBG2_PREFIX_FORMAT + "wtdbg2.log",
    params: args = get_wtdbg2_args_from_wildcards,
            cache_args = lambda wildcards, input, output: "--load-nodes '{cached_nodes}' --load-clips '{cached_clips}' --load-kbm '{cached_kbm}'".format(cached_nodes = input.cached_nodes, cached_clips = input.cached_clips, cached_kbm = input.cached_kbm) if is_wtdbg2_using_cache_from_wildcards(wildcards) else "--dump-kbm '{kbm}'".format(kbm = output.kbm),
            inject_unitigs_args = lambda wildcards, input: "--inject-unitigs '{contigs}'".format(contigs = input.contigs) if is_wtdbg2_injecting_contigs_from_wildcards(wildcards) else "",
            inject_fragment_unitigs_args = lambda wildcards, input: "--inject-fragment-unitigs '{contigs}'".format(contigs = input.fragment_contigs) if is_wtdbg2_injecting_fragment_contigs_from_wildcards(wildcards) else "",
            genome_len_arg = get_genome_len_from_wildcards,
            output_prefix = WTDBG2_PREFIX_FORMAT + "wtdbg2",
            frg_dot_escaped = lambda wildcards, output: ["'" + file + "'" for file in output.frg_dot],
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 100000),
    shell: """
        '{input.binary}' -g {params.genome_len_arg} {params.args} -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' {params.cache_args} {params.inject_unitigs_args} {params.inject_fragment_unitigs_args} 2>&1 | tee '{output.log}'

        if [ ! -z '{input.cached_kbm}' ]; then
            ln -sr '{input.cached_kbm}' '{output.kbm}'
        fi

        if [ ! -z '{input.cached_clips}' ]; then
            ln -sr '{input.cached_clips}' '{output.clips}'
        fi

        for file in {params.frg_dot_escaped}; do
            touch "$file"
        done
    """

rule wtdbg2_consensus:
    input: reads = get_genome_reads_from_wildcards,
           contigs = WTDBG2_PREFIX_FORMAT + "wtdbg2.ctg.lay",
           binary = "external-software/wtdbg2/wtpoa-cns",
    output: consensus = WTDBG2_PREFIX_FORMAT + "wtdbg2.raw.fa",
            completed = touch(WTDBG2_PREFIX_FORMAT + "wtdbg2.raw.fa.completed"),
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 8000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: "{input.binary} -t {threads} -i '{input.contigs}' -fo '{output.consensus}'"

##################
###### Flye ######
##################

localrules: install_flye
rule install_flye:
    output: script = "external-software/Flye/bin/flye",
            directory = directory("external-software/Flye"),
    conda:  "config/conda-install-flye-env.yml"
    shell:  """
        mkdir -p external-software
        cd external-software

        rm -rf Flye
        git clone https://github.com/sebschmi/Flye
        cd Flye
        git checkout 7413f5c39b6c8e9ab9caa564e15e7edd4e727cfd

        make
        """

def get_flye_other_args_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have no assembler arguments: {}".format(arguments))
        cli_arguments = assembler_arguments.get("cli_arguments", None)
        if cli_arguments is None:
            raise Exception("Assembler arguments have no cli arguments: {}".format(assembler_arguments))
        if type(cli_arguments) is not Arguments:
            raise Exception("CLI arguments have type '{}' rather than 'Arguments'.".format(type(cli_arguments)))

        cli_arguments.pop("--pacbio-raw", None)
        cli_arguments.pop("--pacbio-corr", None)
        cli_arguments.pop("--pacbio-hifi", None)
        cli_arguments.pop("--nano-raw", None)
        cli_arguments.pop("--nano-corr", None)

        return cli_arguments.to_argument_string()
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_flye_input_arg_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_arguments = arguments.assembler_arguments()
        if assembler_arguments is None:
            raise Exception("Arguments have not assembler arguments: {}".format(arguments))
        cli_arguments = assembler_arguments.get("cli_arguments", None)
        if cli_arguments is None:
            raise Exception("Assembler arguments have no cli arguments: {}".format(assembler_arguments))

        input_args = []
        if "--pacbio-raw" in cli_arguments:
            input_args.append("--pacbio-raw")
        if "--pacbio-corr" in cli_arguments:
            input_args.append("--pacbio-corr")
        if "--pacbio-hifi" in cli_arguments:
            input_args.append("--pacbio-hifi")
        if "--nano-raw" in cli_arguments:
            input_args.append("--nano-raw")
        if "--nano-corr" in cli_arguments:
            input_args.append("--nano-corr")
        
        if len(input_args) == 0:
            raise Exception("No flye input args given")
        elif len(input_args) == 1:
            return input_args[0]
        elif len(input_args) > 1:
            raise Exception("More than one flye input arg given: {}".format(input_args))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule flye:
    input:  reads = get_genome_reads_from_wildcards,
            script = "external-software/Flye/bin/flye"
    output: contigs = ALGORITHM_PREFIX_FORMAT + "flye/assembly.fasta",
            directory = directory(ALGORITHM_PREFIX_FORMAT + "flye"),
    params: flye_args = get_flye_other_args_from_wildcards,
            flye_input_argument = get_flye_input_arg_from_wildcards,
            genome_len_arg = lambda wildcards: "-g " + get_genome_len_from_wildcards(wildcards),
            output_directory = ALGORITHM_PREFIX_FORMAT + "flye",
            completed = touch(ALGORITHM_PREFIX_FORMAT + "flye/.completed"),
    #conda: "config/conda-flye-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 75000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440, 75000),
               mail_type = "END",
    shell: "'{input.script}' {params.genome_len_arg} {params.flye_args} -t {threads} -o '{params.output_directory}' {params.flye_input_argument} '{input.reads}'"

#####################
###### Hifiasm ######
#####################

rule hifiasm:
    input:  reads = get_genome_reads_from_wildcards,
    output: contigs = os.path.join(HIFIASM_PREFIX_FORMAT, "hifiasm", "assembly.p_ctg.gfa"),
            directory = directory(os.path.join(HIFIASM_PREFIX_FORMAT, "hifiasm")),
    params: output_prefix = os.path.join(HIFIASM_PREFIX_FORMAT, "hifiasm", "assembly"),
    conda:  "config/conda-hifiasm-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 50000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360, 50000),
    shell: "hifiasm -t {threads} -o '{params.output_prefix}' '{input.reads}'"

localrules: hifiasm_gfa_to_fa
rule hifiasm_gfa_to_fa:
    input:  gfa = os.path.join(HIFIASM_PREFIX_FORMAT, "hifiasm", "assembly.p_ctg.gfa"),
    output: fa = os.path.join(HIFIASM_PREFIX_FORMAT, "assembly.p_ctg.fa"),
    run:
            with open(input.gfa, 'r') as input_file, open(output.fa, 'w') as output_file:
                for line in input_file:
                    if line[0] != "S":
                        continue

                    columns = line.split("\t")
                    output_file.write(">{}\n{}\n".format(columns[1], columns[2]))


#############################
###### Read Simulation ######
#############################

def get_perfect_read_simulator_args_from_wildcards(wildcards):
    try:
        read_simulator_arguments = Arguments.from_str(wildcards.read_simulator_arguments)
        cli_arguments = read_simulator_arguments.get("cli_arguments", None)
        if cli_arguments is None:
            return ""

        return cli_arguments.to_argument_string()
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule simulate_perfect_reads:
    input:  reference = GENOME_REFERENCE_FORMAT,
            script = "scripts/simulate_perfect_reads.py",
    output: simulated_reads = safe_format(GENOME_SIMULATED_READS_FORMAT, read_simulator_name = "perfect"),
    log:    log = safe_format(GENOME_SIMULATED_READS_FORMAT, read_simulator_name = "perfect") + ".log",
    params: cli_arguments = get_perfect_read_simulator_args_from_wildcards,
    resources:
        time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
        mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1000),
        queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 1000),
    conda:  "config/conda-simulate-perfect-reads-env.yml",
    shell:  "'{input.script}' {params.cli_arguments} --reference '{input.reference}' --output '{output.simulated_reads}'"

rule simulate_hifi_reads_bbmap:
    input:  reference = GENOME_REFERENCE_FORMAT,
    output: simulated_reads = safe_format(GENOME_SIMULATED_READS_FORMAT, read_simulator_name = "bbmap_hifi"),
    log:    log = safe_format(GENOME_SIMULATED_READS_FORMAT, read_simulator_name = "bbmap_hifi") + ".log",
    params: working_directory = lambda wildcards, output: os.path.dirname(output.simulated_reads),
            mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
    resources:
        time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
        mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
        queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 4000),
    conda:  "config/conda-bbmap-env.yml"
    shell:  """
        REFERENCE=$(realpath -s '{input.reference}')
        OUTPUT=$(realpath -s '{output.simulated_reads}')
        LOG=$(realpath -s '{log.log}')

        cd '{params.working_directory}'
        randomreads.sh build=1 \
        -Xmx{params.mem_mb}m \
        ow=t seed=1 \
        ref="$REFERENCE" \
        simplenames=t \
        pacbio=t pbmin=0.001 pbmax=0.01 \
        coverage=30 paired=f \
        gaussianlength=t \
        minlength=9000 midlength=10000 maxlength=12000 \
        out="$OUTPUT" 2>&1 | tee "$LOG"
        """

#########################################
###### Long Read Input Preparation ######
#########################################

def url_file_format(url):
    parts = url.split('.')
    if parts[-1] == "gz":
        return parts[-2]
    else:
        return parts[-1]

def read_url_file_format(genome):
    try:
        format = url_file_format(genomes[genome]["urls"][0])

        if format == "fasta":
            return "fa"
        elif genomes[genome].get("format", None) == "sra":
            return "sra"
        else:
            return format
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

ruleorder: download_raw_source_reads > download_packed_source_reads > convert_source_reads > extract

localrules: download_raw_source_reads
rule download_raw_source_reads:
    output: file = DATADIR + "downloads/{genome}/reads-{index}/reads.{format}",
            completed = touch(DATADIR + "downloads/{genome}/reads-{index}/reads.{format}.completed"),
    params: url = lambda wildcards: genomes[wildcards.genome]["urls"][int(wildcards.index)],
            url_format = lambda wildcards: read_url_file_format(wildcards.genome),
    wildcard_constraints:
        format = "(fa|bam|sra)",
        index = "\d+",
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output.file}')"

        if [ '{params.url_format}' != '{wildcards.format}' ]; then
            echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
            exit 1
        fi

        wget --progress=dot:mega -O '{output.file}' '{params.url}'
        """

localrules: download_packed_source_reads
rule download_packed_source_reads:
    output: file = DATADIR + "downloads/{genome}/packed-reads-{index}/reads.{format}.gz",
            completed = touch(DATADIR + "downloads/{genome}/packed-reads-{index}/reads.{format}.gz.completed"),
    params: url = lambda wildcards: genomes[wildcards.genome]["urls"][int(wildcards.index)],
            url_format = lambda wildcards: read_url_file_format(wildcards.genome),
            checksum = lambda wildcards: genomes[wildcards.genome]["checksums"][int(wildcards.index)] if "checksums" in genomes[wildcards.genome] else "",
    wildcard_constraints:
        format = "(fa|bam|sra)",
        index = "\d+",
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output.file}')"

        if [ '{params.url_format}' != '{wildcards.format}' ]; then
            echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
            exit 1
        fi

        wget --progress=dot:mega -O '{output.file}' '{params.url}'

        if [ ! -z "{params.checksum}" ]; then
            md5sum -c <<< "{params.checksum} {output.file}"
        else
            echo "Assuming file '{output.file}' was downloaded correctly since no checksum was provided."
        fi
        """

def read_raw_input_file_name(wildcards):
    try:
        genome_name = wildcards.genome
        genome_properties = genomes[genome_name]

        if genome_properties["urls"][0].split('.')[-1] == "gz":
            input_file_name = DATADIR + "downloads/" + genome_name + "/packed-reads-" + wildcards.index + "/reads." + read_url_file_format(wildcards.genome)
        else:
            input_file_name = DATADIR + "downloads/" + genome_name + "/reads-" + wildcards.index + "/reads." + read_url_file_format(wildcards.genome)
        return input_file_name
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule convert_source_reads:
    input: file = read_raw_input_file_name,
    output: file = DATADIR + "downloads/{genome}/reads-{index}/reads.converted.fa",
            completed = touch(DATADIR + "downloads/{genome}/reads-{index}/reads.converted.fa.completed"),
    params: file_format = lambda wildcards: read_url_file_format(wildcards.genome)
    wildcard_constraints:
        index = "\d+",
        genome = "((?!downloads).)*",
    conda: "config/conda-convert-reads-env.yml"
    threads: 1
    shell: """
        if [ '{params.file_format}' == 'bam' ]; then
            samtools fasta '{input.file}' > '{output.file}'
        elif [ '{params.file_format}' == 'sra' ]; then
            fastq-dump --stdout --fasta default '{input.file}' > '{output.file}'
        else
            ln -sr -T '{input.file}' '{output.file}'
        fi
        """

rule combine_reads:
    input: files = lambda wildcards: expand(DATADIR + "downloads/{{genome}}/reads-{index}/reads.converted.fa", index=range(len(genomes[wildcards.genome]["urls"])))
    output: reads = DATADIR + "downloads/{genome}/reads/raw_reads.fa",
            completed = touch(DATADIR + "downloads/{genome}/reads/raw_reads.fa.completed"),
    params: input_list = lambda wildcards, input: "'" + "' '".join(input.files) + "'"
    wildcard_constraints:
        genome = "((?!downloads).)*",
    threads: 1
    resources:
        time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
        queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60),
    shell: "cat {params.input_list} > '{output.reads}'"

rule uniquify_ids:
    input:  reads = DATADIR + "downloads/{genome}/reads/raw_reads.fa",
            script = "scripts/uniquify_fasta_ids.py",
    output: reads = DATADIR + "downloads/{genome}/reads/raw_reads.uniquified.fa",
            log = DATADIR + "downloads/{genome}/uniquify.log",
            completed = touch(DATADIR + "downloads/{genome}/reads/reads.uniquified.fa.completed"),
    wildcard_constraints:
        genome = "((?!downloads).)*",
    conda: "config/conda-uniquify-env.yml"
    threads: 1
    resources:
        time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
        queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60),
    shell: "python3 '{input.script}' '{input.reads}' '{output.reads}' 2>&1 | tee '{output.log}'"

def read_input_file_name(wildcards):
    try:
        genome_name = wildcards.genome
        if genome_name in genomes:
            return DATADIR + "downloads/" + genome_name + "/reads/raw_reads.uniquified.fa"
        elif genome_name in corrected_genomes:
            return DATADIR + "corrected_reads/" + genome_name + "/corrected_reads.fa"
        else:
            sys.exit("genome name not found: " + genome_name)
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: link_reads
rule link_reads:
    input: file = read_input_file_name,
    #input: file = DATADIR + "corrected_reads/{genome}/corrected_reads.fa"
    output: file = GENOME_READS_FORMAT,
    wildcard_constraints:
        genome = "((?!downloads).)*",
    threads: 1
    shell: "ln -sr -T '{input.file}' '{output.file}'"

localrules: download_reference_raw
rule download_reference_raw:
    output: reference = DATADIR + "downloads/{genome}/reference/raw_reference.fa",
            completed = touch(DATADIR + "downloads/{genome}/reference/raw_reference.fa.completed"),
    params: url = lambda wildcards, output: genomes[wildcards.genome]["reference"]
    wildcard_constraints:
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output.reference}')"
        wget --progress=dot:mega -O '{output.reference}' '{params.url}'
        """

localrules: download_reference_gzip
rule download_reference_gzip:
    output: reference = DATADIR + "downloads/{genome}/packed-reference/raw_reference.fa.gz",
            completed = touch(DATADIR + "downloads/{genome}/packed-reference/raw_reference.fa.gz.completed"),
    params: url = lambda wildcards, output: genomes[wildcards.genome]["reference"]
    wildcard_constraints:
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output.reference}')"
        wget --progress=dot:mega -O '{output.reference}' '{params.url}'
        """

def reference_input_file_name(wildcards):
    try:
        genome_name = wildcards.genome
        if corrected_genomes is not None and genome_name in corrected_genomes:
            genome_name = corrected_genomes[genome_name]["source_genome"]
        reference = genomes[genome_name]["reference"]

        if reference.split('.')[-1] == "gz":
            input_file_name = DATADIR + "downloads/" + genome_name + "/packed-reference/raw_reference.fa"
        else:
            input_file_name = DATADIR + "downloads/" + genome_name + "/reference/raw_reference.fa"
        return input_file_name
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: link_reference
rule link_reference:
    input: file = reference_input_file_name,
    output: file = GENOME_REFERENCE_FORMAT,
    wildcard_constraints:
        genome = "((?!downloads).)*",
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output.file}')"
        ln -sr -T '{input.file}' '{output.file}'
        """

#############################
###### Corrected Reads ######
#############################

def correction_read_url_file_format(corrected_genome):
    format = url_file_format(corrected_genomes[corrected_genome]["correction_short_reads"][0])

    if format == "fasta":
        return "fa"
    elif corrected_genomes[corrected_genome].get("format", None) == "sra":
        return "sra"
    else:
        return format

localrules: download_correction_short_reads
rule download_correction_short_reads:
    output: file = DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.{format}",
            completed = touch(DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.{format}.completed"),
    params: url = lambda wildcards: corrected_genomes[wildcards.corrected_genome]["correction_short_reads"][int(wildcards.index)],
            url_format = lambda wildcards: correction_read_url_file_format(wildcards.corrected_genome),
            checksum = lambda wildcards: corrected_genomes[wildcards.corrected_genome]["correction_short_reads_checksum"] if "correction_short_reads_checksum" in corrected_genomes[wildcards.corrected_genome] else "",
    wildcard_constraints:
        format = "(fa|bam|sra)",
        index = "\d+",
        corrected_genome = "((?!corrected_reads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output.file}')"

        if [ '{params.url_format}' != '{wildcards.format}' ]; then
            echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
            exit 1
        fi

        if [ -z "{params.checksum}" ]; then
            wget --progress=dot:mega -O '{output.file}' '{params.url}'
        else
            wget --progress=dot:mega -O '{output.file}' '{params.url}'
            echo "Checksum given, but not supported yet"
            exit 1
        fi
        """

rule convert_correction_short_reads:
    input: file = lambda wildcards: DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads." + correction_read_url_file_format(wildcards.corrected_genome),
    output: file = DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.converted.fa",
            completed = touch(DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.converted.fa.completed"),
    params: file_format = lambda wildcards: correction_read_url_file_format(wildcards.corrected_genome),
    wildcard_constraints:
        index = "\d+",
        genome = "((?!corrected_reads).)*",
    conda: "config/conda-convert-reads-env.yml"
    threads: 1
    resources:
        time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
        queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: """
        if [ '{params.file_format}' == 'bam' ]; then
            samtools fasta '{input.file}' > '{output.file}'
        elif [ '{params.file_format}' == 'sra' ]; then
            fastq-dump --stdout --fasta default '{input.file}' > '{output.file}'
        else
            ln -sr -T '{input.file}' '{output.file}'
        fi
        """

rule combine_correction_short_reads:
    input: files = lambda wildcards: expand(DATADIR + "corrected_reads/{{corrected_genome}}/reads-{index}/reads.converted.fa", index=range(len(corrected_genomes[wildcards.corrected_genome]["correction_short_reads"]))),
    output: reads = CORRECTION_SHORT_READS_FORMAT,
            completed = touch(CORRECTION_SHORT_READS_FORMAT + ".completed"),
    params: input_list = lambda wildcards, input: "'" + "' '".join(input.files) + "'",
    wildcard_constraints:
        genome = "((?!corrected_reads).)*",
    threads: 1
    resources:
        time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
        queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: "cat {params.input_list} > '{output.reads}'"

localrules: install_ratatosk
rule install_ratatosk:
    output: binary = "external-software/Ratatosk/build/src/Ratatosk",
    conda: "config/conda-install-ratatosk-env.yml"
    threads: 1
    shell: """
        mkdir -p external-software
        cd external-software

        rm -r Ratatosk
        git clone --recursive https://github.com/GuillaumeHolley/Ratatosk.git
        cd Ratatosk
        git checkout --recurse-submodules 74ca617afb20a7c24d73d20f2dcdf223db303496

        mkdir build
        cd build
        cmake ..
        make
        """

rule ratatosk:
    input:  correction_short_reads = CORRECTION_SHORT_READS_FORMAT,
            long_reads = lambda wildcards: GENOME_READS_FORMAT.format(genome = corrected_genomes[wildcards.corrected_genome]["source_genome"]),
            binary = "external-software/Ratatosk/build/src/Ratatosk",
    output: corrected_long_reads = DATADIR + "corrected_reads/{corrected_genome}/ratatosk/corrected_reads.fa",
            completed = touch(DATADIR + "corrected_reads/{corrected_genome}/ratatosk/corrected_reads.fa.completed"),
    wildcard_constraints:
        corrected_genome = "((?!/ratatosk).)*",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 80000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 2520),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 2520, 80000),
               mail_type = "END",
    shell: """
        {input.binary} -v -c {threads} -s {input.correction_short_reads} -l {input.long_reads} -o {output.corrected_long_reads}
        """

localrules: select_read_corrector
rule select_read_corrector:
    input:  corrected_reads_from_read_corrector = DATADIR + "corrected_reads/{corrected_genome}/ratatosk/corrected_reads.fa",
    output: corrected_reads = DATADIR + "corrected_reads/{corrected_genome}/corrected_reads.fa",
    wildcard_constraints:
        corrected_genome = "((?!/ratatosk).)*",
    threads: 1
    shell: "ln -sr '{input.corrected_reads_from_read_corrector}' '{output.corrected_reads}'"

###############################
###### Target Generators ######
###############################

def create_experiment_path(experiment):
    return DATADIR + "" + experiment + "/"

def create_report_path(experiment, circularised, k, bcalm2_abundance_min):
    return DATADIR + "" + experiment + "/" + ("circular" if circularised else "linear") + ".k" + str(k) + "-a" + str(bcalm2_abundance_min) + ".report.pdf"

def _generate_read_sim_targets_(experiment, config):
    path = create_experiment_path(experiment)
    for circularised in config["circularised"]:
        yield path + ("circular" if circularised else "linear")

def _generate_bcalm2_parameterisation_targets_(experiment, config):
    for target in _generate_read_sim_targets_(experiment, config):
        for k in config["k"]:
            for bcalm2_abundance_min in config["bcalm2_abundance_min"]:
                    yield target + ".k" + str(k) + "-a" + str(bcalm2_abundance_min)

def _generate_algorithm_targets_(experiment, config):
    for target in _generate_bcalm2_parameterisation_targets_(experiment, config):
        for algorithm in config["algorithm"]:
            yield target + "." + algorithm

def _generate_report_targets_(experiment, config):
        for target in _generate_bcalm2_parameterisation_targets_(experiment, config):
            yield target + ".report.pdf"

def _generate_test_targets_(experiment, config):
        for target in _generate_algorithm_targets_(experiment, config):
            yield target + ".is_tested"

def generate_report_targets():
    for experiment, config in experiments_bcalm2.items():
        for target in _generate_report_targets_(experiment, config):
            yield target

def generate_wtdbg2_report_targets():
    for experiment, config in experiments_wtdbg2.items():
        yield create_experiment_path(experiment) + "wtdbg2.wtdbg2-report.pdf"
    yield DATADIR + "aggregated-report.pdf"

def generate_test_report_targets():
    for experiment, config in tests.items():
        for target in _generate_report_targets_(experiment, config):
            yield target

def generate_test_targets():
    for experiment, config in tests.items():
        for target in _generate_test_targets_(experiment, config):
            yield target

def generate_hamcircuit_targets(amount, n, c):
    yield DATADIR + "hamcircuit/random.0-" + str(amount - 1) + ".n" + str(n) + "-c" + str(c) + ".overallreport"
    for target in generate_hamcircuit_single_report_targets(amount, n, c):
        yield target

def generate_hamcircuit_overall_report_targets(amount, n, c):
    yield "scripts/generate_hamcircuit_overall_report.py"
    for target in generate_hamcircuit_single_report_targets(amount, n, c):
        yield target

def generate_hamcircuit_single_report_targets(amount, n, c):
    for i in range(amount):
        yield DATADIR + "hamcircuit/random-" + str(i) + ".n" + str(n) + "-c" + str(c) + ".report"

######################################
###### Input Genome Preparation ######
######################################

rule separate_linear_and_circular:
    input: filtered = DATADIR + "{genome}/filtered.fna", verified = DATADIR + "{genome}/is_genome_verified.log", binary = PROGRAMDIR + "target/release/cli"
    output: circular = DATADIR + "{genome}/circular.fna", linear = DATADIR + "{genome}/linear.fna", log = DATADIR + "{genome}/separate_linear_and_circular.log"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "cp '{input.filtered}' '{output.linear}'; data/target/release/cli circularise-genome --input '{input.filtered}' 2>&1 --output '{output.circular}' | tee '{output.log}'"

rule verify_genome:
    input: file = DATADIR + "{dir}/filtered.fna", binary = PROGRAMDIR + "target/release/cli"
    output: log = DATADIR + "{dir}/is_genome_verified.log"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: DATADIR + "target/release/cli verify-genome --input '{input.file}' 2>&1 | tee '{output.log}'"

rule filter_genome:
    input: file = DATADIR + "{dir}/raw.fna", binary = PROGRAMDIR + "target/release/cli"
    output: file = DATADIR + "{dir}/filtered.fna", genome_name = DATADIR + "{dir}/name.txt", log = DATADIR + "{dir}/filtered.log"
    params: retain = lambda wildcards: "--retain '" + experiments_bcalm2[wildcards.dir]["filter_retain"] + "'" if "filter_retain" in experiments_bcalm2[wildcards.dir] else ""
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: DATADIR + "target/release/cli filter --input '{input.file}' --output '{output.file}' --extract-name '{output.genome_name}' {params.retain} 2>&1 | tee '{output.log}'"

rule extract:
    input: DATADIR + "{dir}/{file}.gz"
    output: DATADIR + "{dir}/{file}"
    params: working_directory = DATADIR + "{dir}"
    wildcard_constraints:
        file=".*(?<!\.gz)"
        #file=r"^.*([^\.]..|.[^g].|..[^z])$"
    conda: "config/conda-extract-env.yml"
    threads: 1
    resources:
        time_min = 1440,
    shell: "cd '{params.working_directory}'; gunzip -k {wildcards.file}.gz"

#rule extract_dot:
#    input: DATADIR + "{dir}/wtdbg2.3.dot.gz"
#    output: DATADIR + "{dir}/wtdbg2.3.dot"
#    conda: "config/conda-extract-env.yml"
#    shell: "cd 'data/{wildcards.dir}'; gunzip -k wtdbg2.3.dot.gz"

rule download_experiment_file:
    output: DATADIR + "{dir}/raw.fna.gz"
    params: url = lambda wildcards, output: experiments_bcalm2[wildcards.dir]["url"]
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname '{output}')"
        cd 'data/{wildcards.dir}'
        wget --progress=dot:mega -O raw.fna.gz {params.url}
    """

###################
###### QUAST ######
###################

localrules: install_quast
rule install_quast:
    output: script = "external-software/quast/quast.py",
            script_directory = directory("external-software/quast/"),
    threads: 1
    shell: """
    mkdir -p external-software
    cd external-software

    git clone https://github.com/sebschmi/quast
    cd quast
    git checkout 1f2ba4cf40f963f89fbabc9d48e13ebd5c65a77d
    """

rule run_quast_bcalm2:
    input: reads = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa",
        reference = DATADIR + "{dir}/{file}.fna",
        script = "external-software/quast/quast.py",
        script_directory = "external-software/quast/"
    output: report = directory(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.quast")
    conda: "config/conda-quast-env.yml"
    threads: 1
    shell: "{input.script} -t {threads} -o {output.report} -r {input.reference} {input.reads}"

rule run_quast:
    input:  contigs = ALGORITHM_PREFIX_FORMAT + "assembly.fa",
            reference = get_genome_reference_from_wildcards,
            script = "external-software/quast/quast.py",
    output: directory = directory(QUAST_PREFIX_FORMAT),
            eaxmax_csv = QUAST_PREFIX_FORMAT + "aligned_stats/EAxmax_plot.csv",
            completed = touch(QUAST_PREFIX_FORMAT + ".completed"),
    conda: "config/conda-quast-env.yml"
    threads: 4
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 50000),
               cpus = 4,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120),
    shell: "{input.script} -t {threads} --no-html --fragmented --large -o '{output.directory}' -r '{input.reference}' '{input.contigs}'"


##################
###### Rust ######
##################

rule build_rust_release:
    input:  PROGRAMDIR + "is_rust_tested.log",
    output: PROGRAMDIR + "target/release/cli",
    params: programdir = PROGRAMDIR,
    conda: "config/conda-rust-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = 4000,
               cpus = MAX_THREADS,
               time_min = 30,
    shell: "cargo build -j {threads} --release --target-dir '{params.programdir}target' --manifest-path 'implementation/Cargo.toml'"

rule test_rust:
    input:  expand("{source}", source = list(rust_sources))
    output: touch(PROGRAMDIR + "is_rust_tested.log")
    params: programdir = PROGRAMDIR,
    conda: "config/conda-rust-env.yml"
    threads: 4
    resources: mem_mb = 4000,
               cpus = 2,
               time_min = 30,
    shell: "cargo test -j {threads} --target-dir '{params.programdir}target' --manifest-path 'implementation/Cargo.toml' 2>&1 | tee '{output}'"

#####################
###### Testing ######
#####################

rule test:
    input: generate_test_targets()

rule test_algorithm:
    input: DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.is_tested"
    output: touch(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.is_tested")

rule test_single_file:
    input: verify = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.verify",
           deterministic = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.deterministic"
    output: touch(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.is_tested")
    threads: 1
    shell: "cmp --silent {input.verify} {input.deterministic}"

rule make_bcalm_output_deterministic:
    input: file = DATADIR + "{dir}/{file}.bcalm2.fa", script = "scripts/make_bcalm_output_deterministic.py"
    output: file = DATADIR + "{dir}/{file}.bcalm2.fa.deterministic"
    threads: 1
    shell: "python scripts/make_bcalm_output_deterministic.py '{input.file}' '{output.file}'"

rule verify_genome_graph:
    input: file = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = PROGRAMDIR + "target/release/cli"
    output: verification_copy = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.verify", log =  DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.properties", latex = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.graphstatistics"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: PROGRAMDIR + "target/release/cli verify --input '{input.file}' --kmer-size {wildcards.k} --output '{output.verification_copy}' --latex '{output.latex}' 2>&1 | tee '{output.log}.tmp' && mv '{output.log}.tmp' '{output.log}'"

rule selftest:
    conda: "config/conda-selftest-env.yml"
    threads: 1
    shell: "echo \"snakemake $(snakemake --version)\"; conda --version; wget --version"

#######################################
###### Genome Graph Construction ######
#######################################

rule bcalm2:
    input: genome = DATADIR + "{dir}/{file}.fna"
    output: unitigs = DATADIR + "{dir}/{file,(circular|linear)}.k{k,[0-9]+}-a{abundance_min,[0-9]+}.bcalm2.fa",
    #params: tmp = DATADIR + "{dir}/{file,(circular|linear)}.k{k,[0-9]+}-a{abundance_min,[0-9]+}.unitigs.bcalm2-tmp/"
    conda: "config/conda-bcalm2-env.yml"
    threads: workflow.cores
    shell: 
        """
        bcalm -nb-cores {threads} -in '{input.genome}' -out '{output.unitigs}' -kmer-size {wildcards.k} -abundance-min {wildcards.abundance_min}
        mv '{output.unitigs}.unitigs.fa' '{output.unitigs}'
        rm data/{wildcards.dir}/{wildcards.file}.k{wildcards.k}-a{wildcards.abundance_min}.bcalm2.*.glue.*
        """

rule report_bcalm2:
    input: generate_report_targets()

#############################
###### ContigValidator ######
#############################

rule install_sdsl:
    output: dir = "external-software/sdsl-lite"
    conda: "config/conda-contigvalidator-env.yml"
    threads: 1
    shell:
        """
        cd external-software
        git clone https://github.com/simongog/sdsl-lite.git
        cd sdsl-lite
        git checkout v2.1.1
        HOME=`pwd` ./install.sh
        """

rule install_contig_validator:
    input: sdsl = "external-software/sdsl-lite"
    output: dir = directory("external-software/ContigValidator")
    conda: "config/conda-contigvalidator-env.yml"
    threads: 1
    shell: 
        """
        cd external-software
        git clone --recursive https://github.com/mayankpahadia1993/ContigValidator.git
        cd ContigValidator/src
        echo 'count_kmers: count_kmers_kmc' >> Makefile
        sed -i 's\\count_kmers: count_kmers_kmc.cpp KMC/kmc_api/kmc_file.o\\count_kmers_kmc: count_kmers_kmc.cpp KMC/kmc_api/kmc_file.o\\g' Makefile
        LIBRARY_PATH="../../sdsl-lite/lib" CPATH="../../sdsl-lite/include" make
        """

rule run_contig_validator:
    input: cv = "external-software/ContigValidator",
        reads = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa",
        reference = DATADIR + "{dir}/{file}.fna"
    output: result = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.contigvalidator",
        exact_alignments = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.exact"),
        bwa_bam = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.bwa.bam"),
        bwa_bam_bai = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.bwa.bam.bai"),
        fn_kmc_pre = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fn.kmc_pre"),
        fn_kmc_suf = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fn.kmc_suf"),
        fp_kmc_pre = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fp.kmc_pre"),
        fp_kmc_suf = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fp.kmc_suf"),
        kmc_kmc_pre = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.kmc.kmc_pre"),
        kmc_kmc_suf = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.kmc.kmc_suf"),
        tp_kmc_pre = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.tp.kmc_pre"),
        tp_kmc_suf = temp(DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.tp.kmc_suf")
    conda: "config/conda-contigvalidator-env.yml"
    threads: 1
    resources:
      contigvalidator = 1
    shell:
        """
        cd external-software/ContigValidator
        # The abundance-min here has nothing to do with the abundance_min from bcalm2
        bash run.sh -suffixsave 0 -abundance-min 1 -kmer-size {wildcards.k} -r '../../{input.reference}' -a '../../{output.result}' -i '../../{input.reads}'
        """

#####################
###### Bandage ######
#####################

rule download_bcalm2_gfa_converter:
    output: "external-software/scripts/convertToGFA.py"
    conda: "config/conda-download-env.yml"
    threads: 1
    shell:
        """
        mkdir -p external-software/scripts
        cd external-software/scripts
        wget https://raw.githubusercontent.com/GATB/bcalm/v2.2.3/scripts/convertToGFA.py
        chmod u+x convertToGFA.py
        """

rule convert_bcalm2_output_to_gfa:
    input: fa = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa",
        converter = "external-software/scripts/convertToGFA.py"
    output: gfa = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.gfa"
    threads: 1
    shell: "external-software/scripts/convertToGFA.py {input.fa} {output.gfa} {wildcards.k}"

rule bandage:
    input: "{file}.gfa"
    output: "{file}.bandage.png"
    conda: "config/conda-bandage-env.yml"
    threads: 1
    shell: "Bandage image {input} {output} --width 1000 --height 1000"

######################
###### Concorde ######
######################

rule install_concorde:
    output: "external-software/concorde/TSP/concorde"
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
    mkdir -p external-software
    cd external-software

    # Download and unpack concorde
    wget http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz
    tar -xf co031219.tgz
    cd concorde

    # Fix whatever incompatible libraries
    sed -i "s/h->h_addr/h->h_addr_list[0]/g" UTIL/safe_io.c

    # Download QSOpt
    wget https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/PIC/qsopt.PIC.a
    mv qsopt.PIC.a qsopt.a
    wget https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/PIC/qsopt.h

    # Build concorde
    QSOPT=$(pwd)
    CC="gcc" CFLAGS="-g -O3" LDFLAGS="-g -O3" ./configure --with-qsopt="$QSOPT"
    make
    """

########################
###### HamCircuit ######
########################

rule single_hamcircuit_n20_c1_0:
    input: generate_hamcircuit_targets(1, 20, 1.0)
    output: touch(DATADIR + "hamcircuit/tested.1.n20-c1.0.touch")

rule ten_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(10, 20, 1.0)
    output: touch(DATADIR + "hamcircuit/tested.10.n20-c1.0.touch")

rule hundred_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(100, 20, 1.0)
    output: touch(DATADIR + "hamcircuit/tested.100.n20-c1.0.touch")

rule thousand_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(1000, 20, 1.0)
    output: touch(DATADIR + "hamcircuit/tested.1000.n20-c1.0.touch")

rule tenthousand_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(10000, 20, 1.0)
    output: touch(DATADIR + "hamcircuit/tested.10000.n20-c1.0.touch")

rule hundred_hamcircuits_n300_c0_65:
    input: DATADIR + "hamcircuit/tested.100.n300-c0.65.touch"

rule k_hamcircuits_n_c:
    input: lambda wildcards: generate_hamcircuit_targets(int(wildcards.k), int(wildcards.n), float(wildcards.c))
    output: touch(DATADIR + "hamcircuit/tested.{k}.n{n}-c{c}.touch")

rule hundred_hamcircuits_n100_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n100-c{c}.touch", c = [0.6, 0.65, 0.7, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n200_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n200-c{c}.touch", c = [0.6, 0.65, 0.7, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n300_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n300-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n400_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n400-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n500_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n500-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n600_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n600-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_nall_call:
    input: expand(DATADIR + "hamcircuit/tested.100.n{n}-c{c}.touch", n = [100, 200, 300, 400, 500, 600], c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hamcircuit_overall_report:
    input: lambda wildcards: generate_hamcircuit_overall_report_targets(int(wildcards.max) + 1, int(wildcards.n), float(wildcards.c))
    output: DATADIR + "hamcircuit/{name}.0-{max}.n{n}-c{c}.overallreport"
    threads: 1
    shell: "scripts/generate_hamcircuit_overall_report.py 'data/hamcircuit/{wildcards.name}' '{wildcards.max}' '{wildcards.n}' '{wildcards.c}'"

rule hamcircuit_report:
    input: preprocesslog = DATADIR + "hamcircuit/{name}.preprocesslog",
           solution_raw = DATADIR + "hamcircuit/{name}.raw.sol",
           solution_preprocessed = DATADIR + "hamcircuit/{name}.preprocessed.sol",
           tsplog_raw = DATADIR + "hamcircuit/{name}.raw.tsplog",
           tsplog_preprocessed = DATADIR + "hamcircuit/{name}.preprocessed.tsplog",
           script = "scripts/generate_hamcircuit_report.py"
    output: report = DATADIR + "hamcircuit/{name}.report"
    threads: 1
    shell: "'{input.script}' 'data/hamcircuit/{wildcards.name}'"

rule hamcircuit_compute_tsp:
    input: tsp = DATADIR + "hamcircuit/{name}.tsp", binary = "external-software/concorde/TSP/concorde"
    output: solution = DATADIR + "hamcircuit/{name}.sol", tsplog = DATADIR + "hamcircuit/{name}.tsplog"
    threads: workflow.cores
    resources:
      concorde = 1
    shell: """
    cd data/hamcircuit
    cp '{wildcards.name}.tsp' '{wildcards.name}.tsptmp'
    sed -i '0,/EDGE_WEIGHT_SECTION/d' '{wildcards.name}.tsptmp'
    sed -i '/EOF/,$d' '{wildcards.name}.tsptmp'
    LINES=$(wc -l '{wildcards.name}.tsptmp')
    LINES=($LINES)
    LINES=${{LINES[0]}}
    UB=$((LINES * 5 + 1))
    '../../{input.binary}' -u ${{UB}}.5 -o '{wildcards.name}.sol' '{wildcards.name}.tsp' 2>&1 | tee '{wildcards.name}.tsplog'
    """

rule hamcircuit_generate:
    input: binary = PROGRAMDIR + "target/release/cli"
    output: tsp_raw = DATADIR + "hamcircuit/{name}.n{n}-c{c}.raw.tsp", tsp_preprocessed = DATADIR + "hamcircuit/{name}.n{n}-c{c}.preprocessed.tsp", preprocesslog = DATADIR + "hamcircuit/{name}.n{n}-c{c}.preprocesslog"
    shadow: "shallow"
    threads: 1
    shell: """
    '{input.binary}' ham-circuit --input none --random n{wildcards.n}+c{wildcards.c} --output-raw '{output.tsp_raw}.tmp' --output-preprocessed '{output.tsp_preprocessed}.tmp' 2>&1 | tee '{output.preprocesslog}.tmp'
    mv '{output.tsp_raw}.tmp' '{output.tsp_raw}'
    mv '{output.tsp_preprocessed}.tmp' '{output.tsp_preprocessed}'
    mv '{output.preprocesslog}.tmp' '{output.preprocesslog}'
    """

###################################
###### Download requirements ######
###################################

localrules: download_and_prepare
rule download_and_prepare:
    input:  reads = expand(GENOME_READS_FORMAT, genome = genomes.keys()),
            correction_reads = expand(CORRECTION_SHORT_READS_FORMAT, corrected_genome = corrected_genomes.keys()),
            references = expand(GENOME_REFERENCE_FORMAT, genome = genomes.keys()),
            quast = "external-software/quast/quast.py",
            wtdbg2 = "external-software/wtdbg2/wtdbg2",
            rust = PROGRAMDIR + "target/release/cli",
            ratatosk = "external-software/Ratatosk/build/src/Ratatosk",

#rule prepare_wtdbg2:

##############################
###### Download results ######
##############################

rule sync_results:
    conda: "config/conda-rsync-env.yml"
    shell: """
        mkdir -p data/reports
        rsync --verbose --recursive --no-relative --include="*/" --include="report.pdf" --include="aggregated-report.pdf" --exclude="*" turso:'/proj/sebschmi/git/practical-omnitigs/data/reports/' data/reports
        """