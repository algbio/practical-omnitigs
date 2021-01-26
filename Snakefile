import pathlib, itertools, sys

###############################
###### Preprocess Config ######
###############################

print("Preprocessing config", flush = True)

import itertools
import sys
import re

configfile: "config/default.yml"

# Allow to configure to use conda from the config file
if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

workflow.global_resources["contigvalidator"] = 1
workflow.global_resources["concorde"] = 1

MAX_THREADS = 56
print("Setting MAX_THREADS to " + str(MAX_THREADS), flush = True)

# Preprocess experiments configuration
experiments_bcalm2 = config["experiments"]["bcalm2"]
genomes = config["genomes"]
corrected_genomes = config["corrected_genomes"]
reports = config["reports"]
aggregated_reports = config["aggregated_reports"]

class Arguments(dict):
    @staticmethod
    def from_list(arguments):
        result = Arguments()
        for key, value in arguments:
            result[key] = value
        if len(arguments) != len(result):
            sys.exit("Found duplicate arguments: " + str(arguments))
        return result

    @staticmethod
    def from_dict(arguments):
        result = Arguments()
        result.update(arguments)
        return result

    @staticmethod
    def argument_to_str(key, value):
        if value is None or value == True:
            return str(key)
        elif value == False:
            return ""
        else:
            return str(key) + "_" + str(value)

    @staticmethod
    def argument_to_argument_string(key, value):
        if value is None or value == True:
            return str(key)
        elif value == False:
            return ""
        else:
            return str(key) + " " + str(value)

    @staticmethod
    def argument_from_str(string):
        string = string.strip()
        if string == "":
            return ("", False)
        elif "_" in string:
            split = string.split("_")
            return (split[0], split[1])
        else:
            return (string, True)

    def __str__(self):
        result = ""
        once = True
        for argument in self:
            string = Arguments.argument_to_str(argument, self[argument])
            if string != "":
                if once:
                    once = False
                else:
                    result += ":"
                result += string
        return result

    @staticmethod
    def from_str(string):
        string = string.strip()
        result = Arguments()
        if len(string) == 0:
            return result
        for argument in string.split(":"):
            (key, value) = Arguments.argument_from_str(argument)
            result[key] = value
        return result

    def copy(self):
        result = Arguments()
        for key, value in self.items():
            result[key] = value
        return result

    def __key(self):
        return str(self)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Arguments):
            return self.__key() == other.__key()
        return NotImplemented

    def to_argument_string(self):
        result = ""
        once = True
        for key, value in self.items():
            if key.startswith("wtdbg2-") or key.startswith("flye-"):
                continue
            argument = Arguments.argument_to_argument_string(key, value)
            if argument != "":
                if once:
                    once = False
                else:
                    result += " "
                result += argument
        return result

class Algorithm:
    def __init__(self, assembler, arguments):
        self.assembler = assembler
        self.arguments = arguments

    def __str__(self):
        return self.assembler + "::" + str(self.arguments)

    def from_str(string):
        split = string.split("::")
        return Algorithm(split[0], Arguments.from_str(split[1]))

class Experiment:
    def __init__(self, genome, algorithm):
        self.genome = genome
        self.algorithm = algorithm

class Column:
    def __init__(self, shortname, experiment):
        self.shortname = shortname
        self.experiment = experiment

class ReportFile:
    def __init__(self, assembler_name_indices, assembler_argument_map_combination, columns):
        self.columns = columns
        # print("Creating report file with columns: " + str([str(column.experiment.algorithm) for column in columns]))

        self.name = ""
        once = True
        for name, index in assembler_name_indices.items():
            if once:
                once = False
            else:
                self.name += ":::"
            self.name += name + "::" + str(assembler_argument_map_combination[index])

for report_name, report_definition in reports.items():
    report_definition["report_files"] = {}
    assembler_argument_maps = {}

    # Build argument maps for each assembler
    for assembler_name, assembler in report_definition["assemblers"].items():
        for argument_set in assembler["argument_sets"]:
            argument_lists = []
            for key, values in argument_set.items():
                argument_list = []
                for value in values:
                    argument_list.append((key, value))
                argument_lists.append(argument_list)

            argument_maps = []
            for argument_combination in itertools.product(*argument_lists):
                argument_map = Arguments(list(argument_combination))
                argument_maps.append(argument_map)
            assembler_argument_maps.setdefault(assembler_name, set()).update(argument_maps)
        if assembler_name not in assembler_argument_maps or len(assembler_argument_maps[assembler_name]) == 0:
            assembler_argument_maps[assembler_name] = set()
            assembler_argument_maps[assembler_name].add(Arguments())

    # Convert argument maps to list of lists to iterate over their product
    assembler_name_indices = {}
    for index, name in enumerate(assembler_argument_maps.keys()):
        assembler_name_indices[name] = index
    assembler_argument_map_lists = [None] * len(assembler_name_indices)
    for name, index in assembler_name_indices.items():
        assembler_argument_map_list = list(assembler_argument_maps[name])
        assembler_argument_map_lists[index] = assembler_argument_map_list

    # Iterate over the product of argument map lists to create all report files'
    for assembler_argument_map_combination in itertools.product(*assembler_argument_map_lists):
        for genome in report_definition["genomes"]:
            report_file_columns = []
            for column in report_definition["columns"]:
                assembler_name = column["assembler"]
                assembler_arguments = assembler_argument_map_combination[assembler_name_indices[assembler_name]].copy()

                if "arguments" in column:
                    additional_arguments = Arguments.from_dict(column["arguments"])
                    assembler_arguments.update(additional_arguments)
                else:
                    additional_arguments = Arguments()

                algorithm = Algorithm(assembler_name, assembler_arguments)
                experiment = Experiment(genome, algorithm)
                report_file_columns.append(Column(column["shortname"], experiment))
            report_file = ReportFile(assembler_name_indices, assembler_argument_map_combination, report_file_columns)
            report_definition["report_files"][report_file.name] = report_file

import pprint
pp = pprint.PrettyPrinter(indent = 2, width = 200, compact = True)
pp.pprint(reports)
pp.pprint(aggregated_reports)

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

ALGORITHM_PREFIX_FORMAT = "data/{genome}/{algorithm}/"
QUAST_PREFIX_FORMAT = ALGORITHM_PREFIX_FORMAT + "quast/"
REPORT_PREFIX_FORMAT = "data/reports/{report_name}/{report_file_name}"
AGGREGATED_REPORT_PREFIX_FORMAT = "data/reports/{aggregated_report_name}/"

#################################
###### Global report rules ######
#################################

def get_all_report_files():
    result = []
    for report_name, report_definition in reports.items():
        for report_file_name in report_definition["report_files"].keys():
            result.append(REPORT_PREFIX_FORMAT.format(report_name = report_name, report_file_name = report_file_name) + "::::report.pdf")


    for aggregated_report_name in aggregated_reports.keys():
        result.append(AGGREGATED_REPORT_PREFIX_FORMAT.format(aggregated_report_name = aggregated_report_name) + "aggregated-report.pdf")
    return result

rule report_all:
    input: get_all_report_files()
    threads: 1
    resources: mail_type = "END,FAIL,INVALID_DEPEND,REQUEUE"

###############################
###### Report Generation ######
###############################

localrules: create_single_bcalm2_report_tex
rule create_single_bcalm2_report_tex:
    input: genome_name = "data/{dir}/name.txt",
           unitigs = "data/{dir}/{file}.unitigs.tex",
           unitigs_contigvalidator = "data/{dir}/{file}.unitigs.contigvalidator",
           unitigs_quast = "data/{dir}/{file}.unitigs.quast",
           omnitigs = "data/{dir}/{file}.omnitigs.tex",
           omnitigs_contigvalidator = "data/{dir}/{file}.omnitigs.contigvalidator",
           omnitigs_quast = "data/{dir}/{file}.omnitigs.quast",
           trivialomnitigs = "data/{dir}/{file}.trivialomnitigs.tex",
           trivialomnitigs_contigvalidator = "data/{dir}/{file}.trivialomnitigs.contigvalidator",
           trivialomnitigs_quast = "data/{dir}/{file}.trivialomnitigs.quast",
           graphstatistics = "data/{dir}/{file}.bcalm2.graphstatistics",
           bcalm2_bandage = "data/{dir}/{file}.bcalm2.bandage.png",
           script = "scripts/convert_validation_outputs_to_latex.py",
    output: "data/{dir}/{file}.report.tex"
    params: prefix = "data/{dir}/{file}"
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    shell: "python3 scripts/convert_validation_outputs_to_latex.py '{input.genome_name}' '{input.graphstatistics}' '{input.bcalm2_bandage}' 'none' '{output}' uni '{params.prefix}.unitigs' 'Y-to-V' '{params.prefix}.trivialomnitigs' omni '{params.prefix}.omnitigs'"

### Create single report ###

def get_report_file(report_name, report_file_name):
    return reports[report_name]["report_files"][report_file_name]

def get_report_file_quasts(report_name, report_file_name):
    report_file = get_report_file(report_name, report_file_name)
    quasts = []
    for column in report_file.columns:
        genome = column.experiment.genome
        algorithm = column.experiment.algorithm
        quasts.append(QUAST_PREFIX_FORMAT.format(genome = genome, algorithm = algorithm))
    return quasts

def get_report_file_quasts_from_wildcards(wildcards):
    return get_report_file_quasts(wildcards.report_name, wildcards.report_file_name)

def get_report_file_column_shortnames(report_name, report_file_name):
    report_file = get_report_file(report_name, report_file_name)
    shortnames = []
    for column in report_file.columns:
        shortnames.append(column.shortname)
    return shortnames

def get_report_file_column_shortnames_from_wildcards(wildcards):
    return get_report_file_column_shortnames(wildcards.report_name, wildcards.report_file_name)

def get_report_genome_name(report_name, report_file_name):
    report_file = get_report_file(report_name, report_file_name)
    genome_name = None
    for column in report_file.columns:
        genome = column.experiment.genome
        if genome_name is None:
            genome_name = genome
        elif genome_name != genome:
            sys.exit("Report has two separate genomes")

    return genome_name

def get_report_genome_name_from_wildcards(wildcards):
    return get_report_genome_name(wildcards.report_name, wildcards.report_file_name)

def get_single_report_script_column_arguments(report_name, report_file_name):
    report_file = get_report_file(report_name, report_file_name)
    result = ""
    once = True
    for column in report_file.columns:
        if once:
            once = False
        else:
            result += " "

        genome = column.experiment.genome
        algorithm = column.experiment.algorithm
        result += "'" + column.shortname + "' '" + ALGORITHM_PREFIX_FORMAT.format(genome = genome, algorithm = algorithm) + "'"
    return result

def get_single_report_script_column_arguments_from_wildcards(wildcards):
    return get_single_report_script_column_arguments(wildcards.report_name, wildcards.report_file_name)

localrules: create_single_report_tex
rule create_single_report_tex:
    input: quasts = get_report_file_quasts_from_wildcards,
           combined_eaxmax_plot = REPORT_PREFIX_FORMAT + "::::combined_eaxmax_plot.pdf",
           script = "scripts/convert_validation_outputs_to_latex.py",
    output: report = REPORT_PREFIX_FORMAT + "::::report.tex",
    params: genome_name = get_report_genome_name_from_wildcards,
            script_column_arguments = get_single_report_script_column_arguments_from_wildcards,
            name_file = REPORT_PREFIX_FORMAT + "::::name.txt"
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    shell: """
        echo '{wildcards.report_name} {params.genome_name} {wildcards.report_file_name}' > '{params.name_file}'
        python3 '{input.script}' '{params.name_file}' 'none' 'none' '{input.combined_eaxmax_plot}' '{output}' {params.script_column_arguments}
        """

### Create aggregated report ###

def get_aggregated_report_file_maps(aggregated_report_name):
    result = {}
    for report_name in aggregated_reports[aggregated_report_name]["reports"]:
        result.setdefault(report_name, {}).update(reports[report_name]["report_files"])
    return result

def iterate_aggregated_report_file_source_reports(aggregated_report_name):
    aggregated_report_files = get_aggregated_report_file_maps(aggregated_report_name)
    for report_name, report_file_map in aggregated_report_files.items():
        for report_file_name in report_file_map.keys():
            yield (report_name, report_file_name)

def get_aggregated_report_file_source_report_paths(aggregated_report_name):
    result = set()
    for report_name, report_file_name in iterate_aggregated_report_file_source_reports(aggregated_report_name):
        result.add(REPORT_PREFIX_FORMAT.format(report_name = report_name, report_file_name = report_file_name) + "::::report.tex")
    return result

def get_aggregated_report_file_source_report_paths_from_wildcards(wildcards):
    return get_aggregated_report_file_source_report_paths(wildcards.aggregated_report_name)

localrules: create_aggregated_report_tex
rule create_aggregated_report_tex:
    input: source_reports = get_aggregated_report_file_source_report_paths_from_wildcards,
           script = "scripts/create_aggregated_wtdbg2_report.py",
    output: file = AGGREGATED_REPORT_PREFIX_FORMAT + "aggregated-report.tex"
    params: source_reports_arg = lambda wildcards: "' '".join(get_aggregated_report_file_source_report_paths_from_wildcards(wildcards)),
            source_report_names_arg = lambda wildcards: "' '".join([report_name + "/" + report_file_name for report_name, report_file_name in iterate_aggregated_report_file_source_reports(wildcards.aggregated_report_name)]),
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    shell: """
        python3 '{input.script}' --source-reports '{params.source_reports_arg}' --output '{output.file}'
        """

localrules: create_combined_eaxmax_graph
rule create_combined_eaxmax_graph:
    input: quasts = get_report_file_quasts_from_wildcards,
           script = "scripts/create_combined_eaxmax_plot.py",
    output: REPORT_PREFIX_FORMAT + "::::combined_eaxmax_plot.pdf",
    params: input_quasts = lambda wildcards, input: "' '".join([shortname + "' '" + quast for shortname, quast in zip(get_report_file_column_shortnames_from_wildcards(wildcards), input.quasts)])
    conda: "config/conda-seaborn-env.yml"
    threads: 1
    shell: """
        mkdir -p "$(dirname "{output}")"
        python3 '{input.script}' '{params.input_quasts}' '{output}'
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
    input: "data/{subpath}report.tex"
    output: "data/{subpath}report.pdf"
    conda: "config/conda-latex-env.yml"
    threads: 1
    shell: "tectonic {input}"

########################
###### Algorithms ######
########################

### bcalm2 ###

rule compute_omnitigs:
    input: file = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = "data/target/release/cli"
    output: file = "data/{dir}/{file}.k{k}-a{abundance_min}.omnitigs.fa", log = "data/{dir}/{file}.k{k}-a{abundance_min}.omnitigs.fa.log", latex = "data/{dir}/{file}.k{k}-a{abundance_min}.omnitigs.tex"
    threads: 1
    shell: "'{input.binary}' compute-omnitigs --input '{input.file}' --kmer-size {wildcards.k} --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

rule compute_trivial_omnitigs:
    input: file = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = "data/target/release/cli"
    output: file = "data/{dir}/{file}.k{k}-a{abundance_min}.trivialomnitigs.fa", log = "data/{dir}/{file}.k{k}-a{abundance_min}.trivialomnitigs.fa.log", latex = "data/{dir}/{file}.k{k}-a{abundance_min}.trivialomnitigs.tex"
    threads: 1
    shell: "'{input.binary}' compute-trivial-omnitigs --input '{input.file}' --kmer-size {wildcards.k} --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

rule compute_unitigs:
    input: file = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = "data/target/release/cli"
    output: file = "data/{dir}/{file}.k{k}-a{abundance_min}.unitigs.fa", log = "data/{dir}/{file}.k{k}-a{abundance_min}.unitigs.fa.log", latex = "data/{dir}/{file}.k{k}-a{abundance_min}.unitigs.tex"
    threads: 1
    shell: "'{input.binary}' compute-unitigs --input '{input.file}' --kmer-size {wildcards.k} --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

def get_injectable_contigs_rust_cli_command_from_wildcards(wildcards):
    algorithm = Algorithm.from_str(wildcards.algorithm)
    arguments = algorithm.arguments
    if "wtdbg2-inject-unitigs" in arguments:
        return "compute-unitigs"
    elif "wtdbg2-inject-trivial-omnitigs" in arguments:
        return "compute-trivial-omnitigs --non-scc"
    elif "wtdbg2-inject-omnitigs" in arguments:
        return "compute-omnitigs"
    else:
        sys.exit("Missing injection command in wildcards: " + str(wildcards))

rule compute_injectable_contigs_wtdbg2:
    input: nodes = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.3.nodes",
           reads = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.3.reads",
           dot = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.3.dot",
           raw_reads = "data/{genome}/reads.fa",
           binary = "data/target/release/cli"
    output: file = ALGORITHM_PREFIX_FORMAT + "contigwalks",
            log = ALGORITHM_PREFIX_FORMAT + "compute_injectable_contigs.log",
            latex = ALGORITHM_PREFIX_FORMAT + "compute_injectable_contigs.tex"
    params: command = get_injectable_contigs_rust_cli_command_from_wildcards
    threads: 1
    resources: mem_mb = 48000
    shell: "'{input.binary}' {params.command} --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

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
    git checkout 7440d6b5d8ee53a3f5fea9163208b7b295407648
    make
    """

def get_source_genome_properties_from_wildcards(wildcards):
    genome_name = wildcards.genome
    if genome_name in genomes:
        pass
    elif corrected_genomes is not None and genome_name in corrected_genomes:
        genome_name = corrected_genomes[genome_name]["source_genome"]
    else:
        sys.exit("Genome name not found: " + str(genome_name))
    return genomes[genome_name]

def get_assembler_args_from_wildcards(wildcards):
    algorithm = Algorithm.from_str(wildcards.algorithm)
    return algorithm.arguments.to_argument_string()

def get_genome_len_from_wildcards(wildcards):
    genome_properties = get_source_genome_properties_from_wildcards(wildcards)
    return str(genome_properties["genome_length"])

def compute_genome_mem_mb_from_wildcards(wildcards, base_mem_mb):
    genome_properties = get_source_genome_properties_from_wildcards(wildcards)

    if "assembly_mem_factor" in genome_properties:
        return int(float(genome_properties["assembly_mem_factor"]) * float(base_mem_mb))
    else:
        return base_mem_mb

def compute_genome_time_min_from_wildcards(wildcards, base_time_min):
    genome_properties = get_source_genome_properties_from_wildcards(wildcards)

    if "assembly_time_factor" in genome_properties:
        return int(float(genome_properties["assembly_time_factor"]) * float(base_time_min))
    else:
        return base_time_min

def compute_genome_queue_from_wildcards(wildcards, base_time_min):
    time = compute_genome_time_min_from_wildcards(wildcards, base_time_min)
    if time <= 1440:
        return "short"
    elif time <= 1440 * 3:
        return "medium"
    elif time <= 1440 * 14:
        return "long"
    elif time <= 1440 * 60:
        return "extralong"
    else:
        sys.exit("No applicable queue for runtime " + str(time) + " (wildcards: " + str(wildcards) + ")")

rule wtdbg2_complete:
    input: reads = "data/{genome}/reads.fa",
           binary = "external-software/wtdbg2/wtdbg2",
    output: original_nodes = ALGORITHM_PREFIX_FORMAT + "wtdbg2.1.nodes",
            nodes = ALGORITHM_PREFIX_FORMAT + "wtdbg2.3.nodes",
            reads = ALGORITHM_PREFIX_FORMAT + "wtdbg2.3.reads",
            dot = ALGORITHM_PREFIX_FORMAT + "wtdbg2.3.dot.gz",
            clips = ALGORITHM_PREFIX_FORMAT + "wtdbg2.clps",
            kbm = ALGORITHM_PREFIX_FORMAT + "wtdbg2.kbm",
            ctg_lay = ALGORITHM_PREFIX_FORMAT + "wtdbg2.ctg.lay.gz",
            log = ALGORITHM_PREFIX_FORMAT + "wtdbg2.log",
    wildcard_constraints: algorithm = "((?!(--skip-fragment-assembly|wtdbg2-inject-)).)*"
    params: wtdbg2_args = get_assembler_args_from_wildcards,
            genome_len_arg = lambda wildcards: "-g " + get_genome_len_from_wildcards(wildcards),
            output_prefix = ALGORITHM_PREFIX_FORMAT + "wtdbg2",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 60000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: "{input.binary} {params.genome_len_arg} {params.wtdbg2_args} -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --dump-kbm '{output.kbm}' 2>&1 | tee '{output.log}'"

def get_wtdbg2_caching_prefix_from_wildcards(wildcards):
    algorithm = Algorithm.from_str(wildcards.algorithm)
    algorithm.arguments.pop("--skip-fragment-assembly", None)
    algorithm.arguments.pop("wtdbg2-inject-unitigs", None)
    algorithm.arguments.pop("wtdbg2-inject-trivial-omnitigs", None)
    algorithm.arguments.pop("wtdbg2-inject-omnitigs", None)
    algorithm.arguments["wtdbg2-wtdbg2"] = True
    return ALGORITHM_PREFIX_FORMAT.format(genome = "{genome}", algorithm = str(algorithm))

rule wtdbg2_without_fragment_assembly:
    input: reads = "data/{genome}/reads.fa",
           clips = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.clps",
           nodes = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.1.nodes",
           kbm = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.kbm",
           binary = "external-software/wtdbg2/wtdbg2",
    output: original_nodes = ALGORITHM_PREFIX_FORMAT + "wtdbg2.1.nodes",
            nodes = ALGORITHM_PREFIX_FORMAT + "wtdbg2.3.nodes",
            reads = ALGORITHM_PREFIX_FORMAT + "wtdbg2.3.reads",
            dot = ALGORITHM_PREFIX_FORMAT + "wtdbg2.3.dot.gz",
            ctg_lay = ALGORITHM_PREFIX_FORMAT + "wtdbg2.ctg.lay.gz",
            log = ALGORITHM_PREFIX_FORMAT + "wtdbg2.log",
    wildcard_constraints: algorithm = "(.*wtdbg2-wtdbg2.*--skip-fragment-assembly.*|.*--skip-fragment-assembly.*wtdbg2-wtdbg2.*)"
    params: wtdbg2_args = get_assembler_args_from_wildcards,
            genome_len_arg = lambda wildcards: "-g " + get_genome_len_from_wildcards(wildcards),
            output_prefix = ALGORITHM_PREFIX_FORMAT + "wtdbg2",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 48000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: "{input.binary} {params.genome_len_arg} {params.wtdbg2_args} -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --load-nodes '{input.nodes}' --load-clips '{input.clips}' --load-kbm '{input.kbm}' 2>&1 | tee '{output.log}'"

rule wtdbg2_inject_contigs:
    input: reads = "data/{genome}/reads.fa",
           contigs = ALGORITHM_PREFIX_FORMAT + "contigwalks",
           clips = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.clps",
           nodes = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.1.nodes",
           kbm = lambda wildcards: get_wtdbg2_caching_prefix_from_wildcards(wildcards) + "wtdbg2.kbm",
           binary = "external-software/wtdbg2/wtdbg2",
    output: ctg_lay = ALGORITHM_PREFIX_FORMAT + "wtdbg2.ctg.lay.gz",
            log = ALGORITHM_PREFIX_FORMAT + "wtdbg2.log",
    wildcard_constraints: algorithm = ".*wtdbg2-inject-.*"
    params: wtdbg2_args = get_assembler_args_from_wildcards,
            genome_len_arg = lambda wildcards: "-g " + get_genome_len_from_wildcards(wildcards),
            output_prefix = ALGORITHM_PREFIX_FORMAT + "wtdbg2",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 48000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: "{input.binary} {params.genome_len_arg} {params.wtdbg2_args} -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --inject-unitigs '{input.contigs}' --load-nodes '{input.nodes}' --load-clips '{input.clips}' --load-kbm '{input.kbm}' 2>&1 | tee '{output.log}'"

rule wtdbg2_consensus:
    input: reads = "data/{genome}/reads.fa",
           contigs = ALGORITHM_PREFIX_FORMAT + "wtdbg2.ctg.lay",
           binary = "external-software/wtdbg2/wtpoa-cns"
    output: consensus = ALGORITHM_PREFIX_FORMAT + "wtdbg2.raw.fa"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 8000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
    shell: "{input.binary} -t {threads} -i '{input.contigs}' -fo '{output.consensus}'"

##################
###### Flye ######
##################

def get_flye_input_argument_from_wildcards(wildcards):
    algorithm = Algorithm.from_str(wildcards.algorithm)
    if "flye-pacbio-hifi" in algorithm.arguments:
        return "--pacbio-hifi"
    elif "flye-pacbio-corr" in algorithm.arguments:
        return "--pacbio-corr"
    elif "flye-pacbio-raw" in algorithm.arguments:
        return "--pacbio-raw"
    else:
        sys.exit("Missing flye input argument: " + str(algorithm.arguments))

rule flye:
    input: reads = "data/{genome}/reads.fa"
    output: directory = directory(ALGORITHM_PREFIX_FORMAT + "flye/"),
            contigs = ALGORITHM_PREFIX_FORMAT + "flye/assembly.fasta",
    params: flye_args = get_assembler_args_from_wildcards,
            flye_input_argument = get_flye_input_argument_from_wildcards,
            genome_len_arg = lambda wildcards: "-g " + get_genome_len_from_wildcards(wildcards),
    conda: "config/conda-flye-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 250000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440),
               mail_type = "END",
    shell: "flye {params.genome_len_arg} {params.flye_args} -t {threads} -o '{output.directory}' {params.flye_input_argument} '{input.reads}'"


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
    format = url_file_format(genomes[genome]["urls"][0])

    if format == "fasta":
        return "fa"
    elif genomes[genome].get("format", None) == "sra":
        return "sra"
    else:
        return format

ruleorder: download_raw_source_reads > download_packed_source_reads > convert_source_reads > extract

localrules: download_raw_source_reads
rule download_raw_source_reads:
    output: file = "data/downloads/{genome}/reads-{index}.{format}"
    params: url = lambda wildcards: genomes[wildcards.genome]["urls"][int(wildcards.index)],
            url_format = lambda wildcards: read_url_file_format(wildcards.genome)
    wildcard_constraints:
        format = "(fa|bam|sra)",
        index = "\d+",
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p 'data/downloads/{wildcards.genome}'

        if [ '{params.url_format}' != '{wildcards.format}' ]; then
            echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
            exit 1
        fi

        wget --progress=dot:mega -O '{output.file}' '{params.url}'
        """

localrules: download_packed_source_reads
rule download_packed_source_reads:
    output: file = "data/downloads/{genome}/packed-reads-{index}.{format}.gz"
    params: url = lambda wildcards: genomes[wildcards.genome]["urls"][int(wildcards.index)],
            url_format = lambda wildcards: read_url_file_format(wildcards.genome)
    wildcard_constraints:
        format = "(fa|bam|sra)",
        index = "\d+",
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p 'data/downloads/{wildcards.genome}'

        if [ '{params.url_format}' != '{wildcards.format}' ]; then
            echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
            exit 1
        fi

        wget --progress=dot:mega -O '{output.file}' '{params.url}'
        """

def read_raw_input_file_name(wildcards):
    genome_name = wildcards.genome
    genome_properties = genomes[genome_name]

    if genome_properties["urls"][0].split('.')[-1] == "gz":
        input_file_name = "data/downloads/" + genome_name + "/packed-reads-" + wildcards.index + "." + read_url_file_format(wildcards.genome)
    else:
        input_file_name = "data/downloads/" + genome_name + "/reads-" + wildcards.index + "." + read_url_file_format(wildcards.genome)
    return input_file_name

rule convert_source_reads:
    input: file = read_raw_input_file_name
    output: file = "data/downloads/{genome}/reads-{index}.converted.fa"
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
    input: files = lambda wildcards: expand("data/downloads/{{genome}}/reads-{index}.converted.fa", index=range(len(genomes[wildcards.genome]["urls"])))
    output: reads = "data/downloads/{genome}/reads.fa"
    params: input_list = lambda wildcards, input: "'" + "' '".join(input.files) + "'"
    wildcard_constraints:
        genome = "((?!downloads).)*",
    threads: 1
    shell: "cat {params.input_list} > '{output.reads}'"

rule uniquify_ids:
    input: reads = "data/downloads/{genome}/reads.fa", script = "scripts/uniquify_fasta_ids.py"
    output: reads = "data/downloads/{genome}/reads.uniquified.fa", log = "data/downloads/{genome}/uniquify.log"
    wildcard_constraints:
        genome = "((?!downloads).)*",
    conda: "config/conda-uniquify-env.yml"
    threads: 1
    shell: "python3 '{input.script}' '{input.reads}' '{output.reads}' 2>&1 | tee '{output.log}'"

def read_input_file_name(wildcards):
    genome_name = wildcards.genome
    if genome_name in genomes:
        return "data/downloads/" + genome_name + "/reads.uniquified.fa"
    elif genome_name in corrected_genomes:
        return "data/corrected_reads/" + genome_name + "/corrected_reads.fa"
    else:
        sys.exit("genome name not found: " + genome_name)

localrules: link_reads
rule link_reads:
    input: file = read_input_file_name
    #input: file = "data/corrected_reads/{genome}/corrected_reads.fa"
    output: file = "data/{genome}/reads.fa"
    wildcard_constraints:
        genome = "((?!downloads).)*",
    threads: 1
    shell: "ln -sr -T '{input.file}' '{output.file}'"

localrules: download_reference_raw
rule download_reference_raw:
    output: reference = "data/downloads/{genome}/reference.fa"
    params: url = lambda wildcards, output: genomes[wildcards.genome]["reference"]
    wildcard_constraints:
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p 'data/{wildcards.genome}'
        wget --progress=dot:mega -O '{output.reference}' '{params.url}'
        """

localrules: download_reference_gzip
rule download_reference_gzip:
    output: reference = "data/downloads/{genome}/packed-reference.fa.gz"
    params: url = lambda wildcards, output: genomes[wildcards.genome]["reference"]
    wildcard_constraints:
        genome = "((?!downloads).)*",
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p 'data/{wildcards.genome}'
        wget --progress=dot:mega -O '{output.reference}' '{params.url}'
        """

def reference_input_file_name(wildcards):
    genome_name = wildcards.genome
    if corrected_genomes is not None and genome_name in corrected_genomes:
        genome_name = corrected_genomes[genome_name]["source_genome"]
    reference = genomes[genome_name]["reference"]

    if reference.split('.')[-1] == "gz":
        input_file_name = "data/downloads/" + genome_name + "/packed-reference.fa"
    else:
        input_file_name = "data/downloads/" + genome_name + "/reference.fa"
    return input_file_name

localrules: link_reference
rule link_reference:
    input: file = reference_input_file_name
    output: file = "data/{genome}/reference.fa"
    wildcard_constraints:
        genome = "((?!downloads).)*",
    threads: 1
    shell: """
        mkdir -p 'data/{wildcards.genome}'
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
    output: file = "data/corrected_reads/{corrected_genome}/reads-{index}.{format}"
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
        mkdir -p 'data/corrected_reads/{wildcards.corrected_genome}'

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
    input: file = lambda wildcards: "data/corrected_reads/{corrected_genome}/reads-{index}." + correction_read_url_file_format(wildcards.corrected_genome)
    output: file = "data/corrected_reads/{corrected_genome}/reads-{index}.converted.fa"
    params: file_format = lambda wildcards: correction_read_url_file_format(wildcards.corrected_genome)
    wildcard_constraints:
        index = "\d+",
        genome = "((?!corrected_reads).)*",
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

rule combine_correction_short_reads:
    input: files = lambda wildcards: expand("data/corrected_reads/{{corrected_genome}}/reads-{index}.converted.fa", index=range(len(corrected_genomes[wildcards.corrected_genome]["correction_short_reads"])))
    output: reads = "data/corrected_reads/{corrected_genome}/correction_short_reads.fa"
    params: input_list = lambda wildcards, input: "'" + "' '".join(input.files) + "'"
    wildcard_constraints:
        genome = "((?!corrected_reads).)*",
    threads: 1
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
    input: correction_short_reads = "data/corrected_reads/{corrected_genome}/correction_short_reads.fa",
           long_reads = lambda wildcards: "data/downloads/" + corrected_genomes[wildcards.corrected_genome]["source_genome"] + "/reads.uniquified.fa",
           binary = "external-software/Ratatosk/build/src/Ratatosk",
    output: corrected_long_reads = "data/corrected_reads/{corrected_genome}/corrected_reads.fa",
    wildcard_constraints:
        genome = "((?!corrected_reads).)*",
    threads: MAX_THREADS
    resources: mem_mb = 48000,
               cpus = MAX_THREADS,
               time_min = 360,
               mail_type = "END",
    shell: """
        {input.binary} -v -c {threads} -s {input.correction_short_reads} -l {input.long_reads} -o {output.corrected_long_reads}
        """

###############################
###### Target Generators ######
###############################

def create_experiment_path(experiment):
    return "data/" + experiment + "/"

def create_report_path(experiment, circularised, k, bcalm2_abundance_min):
    return "data/" + experiment + "/" + ("circular" if circularised else "linear") + ".k" + str(k) + "-a" + str(bcalm2_abundance_min) + ".report.pdf"

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
    yield "data/aggregated-report.pdf"

def generate_test_report_targets():
    for experiment, config in tests.items():
        for target in _generate_report_targets_(experiment, config):
            yield target

def generate_test_targets():
    for experiment, config in tests.items():
        for target in _generate_test_targets_(experiment, config):
            yield target

def generate_hamcircuit_targets(amount, n, c):
    yield "data/hamcircuit/random.0-" + str(amount - 1) + ".n" + str(n) + "-c" + str(c) + ".overallreport"
    for target in generate_hamcircuit_single_report_targets(amount, n, c):
        yield target

def generate_hamcircuit_overall_report_targets(amount, n, c):
    yield "scripts/generate_hamcircuit_overall_report.py"
    for target in generate_hamcircuit_single_report_targets(amount, n, c):
        yield target

def generate_hamcircuit_single_report_targets(amount, n, c):
    for i in range(amount):
        yield "data/hamcircuit/random-" + str(i) + ".n" + str(n) + "-c" + str(c) + ".report"

######################################
###### Input Genome Preparation ######
######################################

rule separate_linear_and_circular:
    input: filtered = "data/{genome}/filtered.fna", verified = "data/{genome}/is_genome_verified.log", binary = "data/target/release/cli"
    output: circular = "data/{genome}/circular.fna", linear = "data/{genome}/linear.fna", log = "data/{genome}/separate_linear_and_circular.log"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "cp '{input.filtered}' '{output.linear}'; data/target/release/cli circularise-genome --input '{input.filtered}' 2>&1 --output '{output.circular}' | tee '{output.log}'"

rule verify_genome:
    input: file = "data/{dir}/filtered.fna", binary = "data/target/release/cli"
    output: log = "data/{dir}/is_genome_verified.log"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "data/target/release/cli verify-genome --input '{input.file}' 2>&1 | tee '{output.log}'"

rule filter_genome:
    input: file = "data/{dir}/raw.fna", binary = "data/target/release/cli"
    output: file = "data/{dir}/filtered.fna", genome_name = "data/{dir}/name.txt", log = "data/{dir}/filtered.log"
    params: retain = lambda wildcards: "--retain '" + experiments_bcalm2[wildcards.dir]["filter_retain"] + "'" if "filter_retain" in experiments_bcalm2[wildcards.dir] else ""
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "data/target/release/cli filter --input '{input.file}' --output '{output.file}' --extract-name '{output.genome_name}' {params.retain} 2>&1 | tee '{output.log}'"

rule extract:
    input: "data/{dir}/{file}.gz"
    output: "data/{dir}/{file}"
    wildcard_constraints:
        file=".*(?<!\.gz)"
        #file=r"^.*([^\.]..|.[^g].|..[^z])$"
    conda: "config/conda-extract-env.yml"
    threads: 1
    shell: "cd 'data/{wildcards.dir}'; gunzip -k {wildcards.file}.gz"

#rule extract_dot:
#    input: "data/{dir}/wtdbg2.3.dot.gz"
#    output: "data/{dir}/wtdbg2.3.dot"
#    conda: "config/conda-extract-env.yml"
#    shell: "cd 'data/{wildcards.dir}'; gunzip -k wtdbg2.3.dot.gz"

rule download_experiment_file:
    output: "data/{dir}/raw.fna.gz"
    params: url = lambda wildcards, output: experiments_bcalm2[wildcards.dir]["url"]
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: "mkdir -p 'data/{wildcards.dir}'; cd 'data/{wildcards.dir}'; wget --progress=dot:mega -O raw.fna.gz {params.url}"

###################
###### QUAST ######
###################

localrules: install_quast
rule install_quast:
    output: script = "external-software/quast/quast.py", script_directory = directory("external-software/quast/")
    conda: "config/conda-install-env.yml"
    threads: 1
    shell: """
    mkdir -p external-software
    cd external-software

    git clone https://github.com/sebschmi/quast
    cd quast
    git checkout 673a78601ac2453b2d994c20f6d998a05ca88fa9
    """

rule run_quast_bcalm2:
    input: reads = "data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa",
        reference = "data/{dir}/{file}.fna",
        script = "external-software/quast/quast.py",
        script_directory = "external-software/quast/"
    output: report = directory("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.quast")
    conda: "config/conda-quast-env.yml"
    threads: 1
    shell: "{input.script} -t {threads} -o {output.report} -r {input.reference} {input.reads}"

rule run_quast_wtdbg2:
    input: contigs = ALGORITHM_PREFIX_FORMAT + "wtdbg2.raw.fa",
        reference = "data/{genome}/reference.fa",
        script = "external-software/quast/quast.py"
    output: report = directory(QUAST_PREFIX_FORMAT)
    wildcard_constraints: algorithm = "wtdbg2::.*"
    conda: "config/conda-quast-env.yml"
    threads: 8
    resources: mem_mb = 12000,
               cpus = 4,
               time_min = 60,
    shell: "{input.script} -t {threads} --fragmented --large -o {output.report} -r {input.reference} {input.contigs}"

rule run_quast_flye:
    input: contigs = ALGORITHM_PREFIX_FORMAT + "flye/assembly.fasta",
        reference = "data/{genome}/reference.fa",
        script = "external-software/quast/quast.py",
        script_directory = "external-software/quast/"
    output: report = directory(QUAST_PREFIX_FORMAT)
    wildcard_constraints: algorithm = "flye::.*"
    conda: "config/conda-quast-env.yml"
    threads: 8
    resources: mem_mb = 12000,
               cpus = 4,
               time_min = 60,
    shell: "{input.script} -t {threads} --fragmented --large -o {output.report} -r {input.reference} {input.contigs}"


##################
###### Rust ######
##################

rule build_rust_release:
    input: "data/is_rust_tested.log"
    output: "data/target/release/cli"
    conda: "config/conda-rust-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = 4000,
               cpus = MAX_THREADS,
               time_min = 30,
    shell: "cargo build -j {threads} --release --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

rule test_rust:
    input: expand("{source}", source = list(rust_sources))
    output: touch("data/is_rust_tested.log")
    conda: "config/conda-rust-env.yml"
    threads: 4
    resources: mem_mb = 4000,
               cpus = 2,
               time_min = 30,
    shell: "cargo test -j {threads} --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml' 2>&1 | tee '{output}'"

#####################
###### Testing ######
#####################

rule test:
    input: generate_test_targets()

rule test_algorithm:
    input: "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.is_tested"
    output: touch("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.is_tested")

rule test_single_file:
    input: verify = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.verify",
           deterministic = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.deterministic"
    output: touch("data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.is_tested")
    threads: 1
    shell: "cmp --silent {input.verify} {input.deterministic}"

rule make_bcalm_output_deterministic:
    input: file = "data/{dir}/{file}.bcalm2.fa", script = "scripts/make_bcalm_output_deterministic.py"
    output: file = "data/{dir}/{file}.bcalm2.fa.deterministic"
    threads: 1
    shell: "python scripts/make_bcalm_output_deterministic.py '{input.file}' '{output.file}'"

rule verify_genome_graph:
    input: file = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa", binary = "data/target/release/cli"
    output: verification_copy = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.verify", log =  "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa.properties", latex = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.graphstatistics"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "data/target/release/cli verify --input '{input.file}' --kmer-size {wildcards.k} --output '{output.verification_copy}' --latex '{output.latex}' 2>&1 | tee '{output.log}.tmp' && mv '{output.log}.tmp' '{output.log}'"

rule selftest:
    conda: "config/conda-selftest-env.yml"
    threads: 1
    shell: "echo \"snakemake $(snakemake --version)\"; conda --version; wget --version"

#######################################
###### Genome Graph Construction ######
#######################################

rule bcalm2:
    input: genome = "data/{dir}/{file}.fna"
    output: unitigs = "data/{dir}/{file,(circular|linear)}.k{k,[0-9]+}-a{abundance_min,[0-9]+}.bcalm2.fa",
    #params: tmp = "data/{dir}/{file,(circular|linear)}.k{k,[0-9]+}-a{abundance_min,[0-9]+}.unitigs.bcalm2-tmp/"
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
        reads = "data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa",
        reference = "data/{dir}/{file}.fna"
    output: result = "data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.contigvalidator",
        exact_alignments = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.exact"),
        bwa_bam = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.bwa.bam"),
        bwa_bam_bai = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.bwa.bam.bai"),
        fn_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fn.kmc_pre"),
        fn_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fn.kmc_suf"),
        fp_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fp.kmc_pre"),
        fp_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.fp.kmc_suf"),
        kmc_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.kmc.kmc_pre"),
        kmc_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.kmc.kmc_suf"),
        tp_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.tp.kmc_pre"),
        tp_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa.tp.kmc_suf")
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
    input: fa = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa",
        converter = "external-software/scripts/convertToGFA.py"
    output: gfa = "data/{dir}/{file}.k{k}-a{abundance_min}.bcalm2.gfa"
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
    output: touch("data/hamcircuit/tested.1.n20-c1.0.touch")

rule ten_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(10, 20, 1.0)
    output: touch("data/hamcircuit/tested.10.n20-c1.0.touch")

rule hundred_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(100, 20, 1.0)
    output: touch("data/hamcircuit/tested.100.n20-c1.0.touch")

rule thousand_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(1000, 20, 1.0)
    output: touch("data/hamcircuit/tested.1000.n20-c1.0.touch")

rule tenthousand_hamcircuits_n20_c1_0:
    input: generate_hamcircuit_targets(10000, 20, 1.0)
    output: touch("data/hamcircuit/tested.10000.n20-c1.0.touch")

rule hundred_hamcircuits_n300_c0_65:
    input: "data/hamcircuit/tested.100.n300-c0.65.touch"

rule k_hamcircuits_n_c:
    input: lambda wildcards: generate_hamcircuit_targets(int(wildcards.k), int(wildcards.n), float(wildcards.c))
    output: touch("data/hamcircuit/tested.{k}.n{n}-c{c}.touch")

rule hundred_hamcircuits_n100_call:
    input: expand("data/hamcircuit/tested.100.n100-c{c}.touch", c = [0.6, 0.65, 0.7, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n200_call:
    input: expand("data/hamcircuit/tested.100.n200-c{c}.touch", c = [0.6, 0.65, 0.7, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n300_call:
    input: expand("data/hamcircuit/tested.100.n300-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n400_call:
    input: expand("data/hamcircuit/tested.100.n400-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n500_call:
    input: expand("data/hamcircuit/tested.100.n500-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_n600_call:
    input: expand("data/hamcircuit/tested.100.n600-c{c}.touch", c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hundred_hamcircuits_nall_call:
    input: expand("data/hamcircuit/tested.100.n{n}-c{c}.touch", n = [100, 200, 300, 400, 500, 600], c = [0.65, 0.7, 0.75, 0.8, 0.9, 1.0])

rule hamcircuit_overall_report:
    input: lambda wildcards: generate_hamcircuit_overall_report_targets(int(wildcards.max) + 1, int(wildcards.n), float(wildcards.c))
    output: "data/hamcircuit/{name}.0-{max}.n{n}-c{c}.overallreport"
    threads: 1
    shell: "scripts/generate_hamcircuit_overall_report.py 'data/hamcircuit/{wildcards.name}' '{wildcards.max}' '{wildcards.n}' '{wildcards.c}'"

rule hamcircuit_report:
    input: preprocesslog = "data/hamcircuit/{name}.preprocesslog",
           solution_raw = "data/hamcircuit/{name}.raw.sol",
           solution_preprocessed = "data/hamcircuit/{name}.preprocessed.sol",
           tsplog_raw = "data/hamcircuit/{name}.raw.tsplog",
           tsplog_preprocessed = "data/hamcircuit/{name}.preprocessed.tsplog",
           script = "scripts/generate_hamcircuit_report.py"
    output: report = "data/hamcircuit/{name}.report"
    threads: 1
    shell: "'{input.script}' 'data/hamcircuit/{wildcards.name}'"

rule hamcircuit_compute_tsp:
    input: tsp = "data/hamcircuit/{name}.tsp", binary = "external-software/concorde/TSP/concorde"
    output: solution = "data/hamcircuit/{name}.sol", tsplog = "data/hamcircuit/{name}.tsplog"
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
    input: binary = "data/target/release/cli"
    output: tsp_raw = "data/hamcircuit/{name}.n{n}-c{c}.raw.tsp", tsp_preprocessed = "data/hamcircuit/{name}.n{n}-c{c}.preprocessed.tsp", preprocesslog = "data/hamcircuit/{name}.n{n}-c{c}.preprocesslog"
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

localrules: download_all
rule download_all:
    input: reads = [file for genome in genomes.keys() for file in expand("data/downloads/{genome}/reads-{index}.converted.fa", genome=[genome], index=range(len(genomes[genome]["urls"])))],
           correction_short_reads = [file for corrected_genome in corrected_genomes.keys() for file in
                                     expand("data/corrected_reads/{corrected_genome}/reads-{index}.{file_type}",
                                            corrected_genome=[corrected_genome],
                                            index=range(len(corrected_genomes[corrected_genome]["correction_short_reads"])),
                                            file_type=correction_read_url_file_format(corrected_genome))] if corrected_genomes is not None else [],
           references = expand("data/{genome}/reference.fa", genome=genomes.keys()),
           quast = "external-software/quast/quast.py",
           wtdbg2 = "external-software/wtdbg2/wtdbg2",
           rust = "data/target/release/cli",
           ratatosk = "external-software/Ratatosk/build/src/Ratatosk",

#rule prepare_wtdbg2:

##############################
###### Download results ######
##############################

rule sync_results:
    conda: "config/conda-rsync-env.yml"
    shell: """
        mkdir -p data/reports
        rsync --verbose --recursive --relative turso:'/proj/sebschmi/git/practical-omnitigs/data/reports/./*/*report.pdf' data/reports
        """