import itertools
import sys
import re
import json
import traceback
import pathlib
import parse
from urllib.parse import urlparse
import datetime

###############################
###### Preprocess Config ######
###############################

print("Preprocessing config", flush = True)

configfile: "config/default.yml"

# Allow to configure to use conda from the config file
if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

workflow.global_resources["contigvalidator"] = 1

DATADIR = "data"
if "datadir" in config:
    DATADIR = config["datadir"]

print("DATADIR: {}".format(DATADIR))

MAX_THREADS = 56
print(f"Setting MAX_THREADS to {MAX_THREADS}", flush = True)
BUILD_THREADS = 7
print(f"Setting BUILD_THREADS to {BUILD_THREADS}", flush = True)

# Preprocess experiments configuration
genomes = config["genomes"]
corrected_genomes = config["corrected_genomes"]
reports = config["reports"]
aggregated_reports = config["aggregated_reports"]
if aggregated_reports is None:
    aggregated_reports = {}

from report_file_parser import *

DEFAULT_ASSEMBLER_ARGUMENTS = {
    "wtdbg2": {
        "skip_fragment_assembly": "no",
        "fragment_correction_steps": "all",
        "tig_injection": "none",
        "frg_injection": "none",
        "frg_injection_stage": "none",
        "hodeco_consensus": "none",
    },
    "mdbg": {
        "mdbg_mode": "multik",
    },
    "hifiasm": {
        "contig_algorithm": "builtin",
    },
}

for report_name, report_definition in reports.items():
    argument_matrix = ArgumentMatrix(report_definition.setdefault("argument_matrix", {}))
    report_definition["argument_matrix"] = argument_matrix # use correct type

    # print("Matrix of {} has length {}".format(report_name, len(argument_matrix)))
    # entries = list(iter(argument_matrix))
    # print("Entries in matrix:")
    # for entry in entries:
    #     print(entry)

    for arguments in argument_matrix:
        if "genome" in arguments:
            genome_arguments = Arguments.from_dict(genomes[arguments["genome"]].setdefault("genome_arguments", {}))
            genome_arguments.update(arguments)
            arguments = genome_arguments # update method works in-place

        arguments.setdefault("read_source", "real")
        arguments.setdefault("read_simulation_model_source", "none")
        arguments.setdefault("read_downsampling_factor", "none")
        arguments.setdefault("homopolymer_compression", "none")
        arguments.setdefault("uniquify_ids", "no")
        arguments.setdefault("assembler", None)
        arguments.setdefault("assembler_arguments", None)
        arguments.setdefault("quast_mode", "normal")
        arguments.setdefault("filter_nw", "no")
        arguments.setdefault("retain_cm", "no")
        arguments.setdefault("filter_plasmids", "no")

        # hisim produces reads with duplicate ids, so we always need to uniquify those
        if arguments["read_source"].startswith("hisim_"):
            arguments["uniquify_ids"] = "yes"

        columns = []
        for column_definition in report_definition["columns"]:
            column_arguments = Arguments.from_dict(column_definition)
            column = Column(arguments, column_arguments)

            assembler_name = column.arguments.assembler_name()
            if assembler_name in DEFAULT_ASSEMBLER_ARGUMENTS:
                #print(f"assembler_name: {assembler_name}")
                assembler_arguments = column.arguments.assembler_arguments()
                #print(f"assembler_arguments: {assembler_arguments}")
                for key, value in DEFAULT_ASSEMBLER_ARGUMENTS[assembler_name].items():
                    assembler_arguments.setdefault(key, value)
                #print(f"assembler_arguments: {assembler_arguments}")

            #print(column, flush = True)
            columns.append(column)

        report_file = ReportFile(arguments, columns)
        report_definition.setdefault("report_files", {})[report_file.name] = report_file

#import pprint
#pp = pprint.PrettyPrinter(indent = 2, width = 200, compact = True)
#pp.pprint(reports)
#pp.pprint(aggregated_reports)

#print("Exiting for debugging")
#sys.exit(0)

# Collect all rust sources
RUST_SOURCES = list(map(str, itertools.chain(pathlib.Path('implementation').glob('**/Cargo.toml'), pathlib.Path('implementation').glob('**/*.rs'))))

import datetime
today = datetime.date.today().isoformat()

print("Finished config preprocessing", flush = True)

#########################
###### Directories ######
#########################

EXTERNAL_SOFTWARE_ROOTDIR = os.path.join(DATADIR, "external-software")
DOWNLOAD_ROOTDIR = os.path.join(DATADIR, "downloads")
GENOME_ROOTDIR = os.path.join(DATADIR, "genomes")
SIMULATION_ROOTDIR = os.path.join(DATADIR, "simulation")
ASSEMBLY_ROOTDIR = os.path.join(DATADIR, "assembly")
EVALUATION_ROOTDIR = os.path.join(DATADIR, "evaluation")
REPORT_ROOTDIR = os.path.join(DATADIR, "reports")

# Genomes

GENOME_ARGUMENT_STRING = "g{genome}-h{homopolymer_compression}"
GENOME_REFERENCE_ARGUMENT_STRING = "reference-n{filter_nw}-c{retain_cm}-p{filter_plasmids}"
GENOME_REFERENCE_SUBDIR = os.path.join(GENOME_ARGUMENT_STRING, GENOME_REFERENCE_ARGUMENT_STRING)
GENOME_REFERENCE = os.path.join(GENOME_ROOTDIR, GENOME_REFERENCE_SUBDIR, "reference.fa")
UNFILTERED_GENOME_REFERENCE = safe_format(GENOME_REFERENCE, filter_nw = "no", retain_cm = "no", filter_plasmids = "no")
GENOME_REFERENCE_LENGTH = os.path.join(GENOME_ROOTDIR, GENOME_REFERENCE_SUBDIR, "reference_length.txt")
UNFILTERED_GENOME_REFERENCE_LENGTH = safe_format(GENOME_REFERENCE_LENGTH, filter_nw = "no", retain_cm = "no", filter_plasmids = "no")
GENOME_READS_ARGUMENT_STRING = "reads-s{read_source}-m{read_simulation_model_source}-r{read_downsampling_factor}-u{uniquify_ids}"
GENOME_READS_SUBDIR = os.path.join(GENOME_ARGUMENT_STRING, GENOME_READS_ARGUMENT_STRING)
GENOME_READS_DIR = os.path.join(GENOME_ROOTDIR, GENOME_READS_SUBDIR)
GENOME_READS = os.path.join(GENOME_READS_DIR, "reads.fa")
GENOME_SINGLE_LINE_READS = os.path.join(GENOME_ROOTDIR, GENOME_READS_SUBDIR, "single_line_reads.fa")
UNIQUIFY_IDS_LOG = os.path.join(GENOME_ROOTDIR, GENOME_READS_SUBDIR, "uniquify_ids.log")

HOCO_READS_LOG = os.path.join(GENOME_ROOTDIR, GENOME_READS_SUBDIR, "hoco.log")
HOCO_REFERENCE_LOG = os.path.join(GENOME_ROOTDIR, GENOME_REFERENCE_SUBDIR, "hoco.log")

# Simulations

FASTK_ARGUMENT_STRING = "k{fastk_k}"
FASTK_OUTPUT_DIR = os.path.join(SIMULATION_ROOTDIR, "fastk", GENOME_READS_SUBDIR, FASTK_ARGUMENT_STRING)
FASTK_INPUT_READS = os.path.join(FASTK_OUTPUT_DIR, "reads.fa")
FASTK_TABLE = os.path.join(FASTK_OUTPUT_DIR, "reads.ktab")
FASTK_PROFILE = os.path.join(FASTK_OUTPUT_DIR, "reads.prof")
FASTK_HIST = os.path.join(FASTK_OUTPUT_DIR, "reads.hist")
FASTK_SYMMETRIC_TABLE = os.path.join(FASTK_OUTPUT_DIR, "symmetric_reads", "symmetric_reads.ktab")
FASTK_HISTEX_EVALUATION = os.path.join(EVALUATION_ROOTDIR, "histex", GENOME_READS_SUBDIR, FASTK_ARGUMENT_STRING, "histogram.txt")

HISIM_ARGUMENT_STRING = "a{himodel_kmer_threshold}-l{himodel_min_valid}-h{himodel_max_valid}"
HIMODEL_OUTPUT_DIR = os.path.join(SIMULATION_ROOTDIR, "himodel", GENOME_READS_SUBDIR, FASTK_ARGUMENT_STRING, HISIM_ARGUMENT_STRING)
HIMODEL_INPUT_SYMMETRIC_TABLE = os.path.join(HIMODEL_OUTPUT_DIR, "reads.ktab")
HIMODEL_INPUT_PROFILE = os.path.join(HIMODEL_OUTPUT_DIR, "reads.prof")
HIMODEL_MODEL = os.path.join(HIMODEL_OUTPUT_DIR, "reads.model")
HISIM_HAPLOTYPE = os.path.join(GENOME_READS_DIR, "reads.hap{haplotype_index}.fasta")

# Assemblies

ASSEMBLY_ARGUMENT_STRING = "a{assembler}--{assembler_arguments}-"
ASSEMBLY_SUBDIR = os.path.join(GENOME_READS_SUBDIR, ASSEMBLY_ARGUMENT_STRING)
ASSEMBLY_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR)
ASSEMBLED_CONTIGS = os.path.join(ASSEMBLY_OUTPUT_DIR, "contigs.fa")
ASSEMBLY_LOG = os.path.join(ASSEMBLY_OUTPUT_DIR, "assembly.log")
ASSEMBLER_ARGUMENT_STRINGS = {}

WTDBG2_ARGUMENT_STRING = "m{wtdbg2_mode}-s{skip_fragment_assembly}-f{fragment_correction_steps}-t{tig_injection}-f{frg_injection}-s{frg_injection_stage}-d{hodeco_consensus}-c{retain_cm}"
ASSEMBLER_ARGUMENT_STRINGS["wtdbg2"] = WTDBG2_ARGUMENT_STRING
WTDBG2_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "wtdbg2", assembler_arguments = WTDBG2_ARGUMENT_STRING)
WTDBG2_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, WTDBG2_SUBDIR)
WTDBG2_INJECTABLE_CONTIG_DIR = os.path.join(WTDBG2_OUTPUT_DIR, "injectable_contigs")
WTDBG2_INJECTABLE_FRAGMENT_CONTIG_DIR = os.path.join(WTDBG2_OUTPUT_DIR, "injectable_fragment_contigs")
WTDBG2_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "wtdbg2")
WTDBG2_LOG = os.path.join(WTDBG2_OUTPUT_DIR, "assembly.log")
WTDBG2_EXTRACT_LOG = os.path.join(WTDBG2_OUTPUT_DIR, "extract.{subfile}.log")
WTDBG2_CONSENSUS_LOG = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2_consensus.log")
WTDBG2_ASSEMBLED_CONTIGS = os.path.join(WTDBG2_OUTPUT_DIR, "contigs.fa")

FLYE_ARGUMENT_STRING = "m{flye_mode}-c{retain_cm}"
ASSEMBLER_ARGUMENT_STRINGS["flye"] = FLYE_ARGUMENT_STRING
FLYE_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "flye", assembler_arguments = FLYE_ARGUMENT_STRING)
FLYE_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, FLYE_SUBDIR)
FLYE_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "flye")
FLYE_LOG = os.path.join(FLYE_OUTPUT_DIR, "assembly.log")
FLYE_ASSEMBLED_CONTIGS = os.path.join(FLYE_OUTPUT_DIR, "contigs.fa")

HIFIASM_ARGUMENT_STRING = "c{contig_algorithm}"
ASSEMBLER_ARGUMENT_STRINGS["hifiasm"] = HIFIASM_ARGUMENT_STRING
HIFIASM_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "hifiasm", assembler_arguments = HIFIASM_ARGUMENT_STRING)
HIFIASM_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, HIFIASM_SUBDIR)
HIFIASM_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "hifiasm")
HIFIASM_LOG = os.path.join(HIFIASM_OUTPUT_DIR, "assembly.log")
HIFIASM_ASSEMBLED_CONTIGS = os.path.join(HIFIASM_OUTPUT_DIR, "contigs.fa")

MDBG_ARGUMENT_STRING = "m{mdbg_mode}"
ASSEMBLER_ARGUMENT_STRINGS["mdbg"] = MDBG_ARGUMENT_STRING
MDBG_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "mdbg", assembler_arguments = MDBG_ARGUMENT_STRING)
MDBG_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, MDBG_SUBDIR)
MDBG_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "mdbg")
MDBG_LOG = os.path.join(MDBG_OUTPUT_DIR, "assembly.log")
MDBG_ASSEMBLED_CONTIGS = os.path.join(MDBG_OUTPUT_DIR, "contigs.fa")

LJA_ARGUMENT_STRING = "none"
ASSEMBLER_ARGUMENT_STRINGS["lja"] = LJA_ARGUMENT_STRING
LJA_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "lja", assembler_arguments = LJA_ARGUMENT_STRING)
LJA_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, LJA_SUBDIR)
LJA_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "lja")
LJA_LOG = os.path.join(LJA_OUTPUT_DIR, "assembly.log")
LJA_ASSEMBLED_CONTIGS = os.path.join(LJA_OUTPUT_DIR, "contigs.fa")

CANU_ARGUMENT_STRING = "c{retain_cm}"
ASSEMBLER_ARGUMENT_STRINGS["canu"] = CANU_ARGUMENT_STRING
CANU_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "canu", assembler_arguments = CANU_ARGUMENT_STRING)
CANU_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, CANU_SUBDIR)
CANU_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "canu")
CANU_LOG = os.path.join(CANU_OUTPUT_DIR, "assembly.log")
CANU_ASSEMBLED_CONTIGS = os.path.join(CANU_OUTPUT_DIR, "contigs.fa")

REFASM_ARGUMENT_STRING = "none"
ASSEMBLER_ARGUMENT_STRINGS["refasm"] = REFASM_ARGUMENT_STRING
REFASM_SUBDIR = safe_format(ASSEMBLY_SUBDIR, assembler = "refasm", assembler_arguments = REFASM_ARGUMENT_STRING)
REFASM_OUTPUT_DIR = os.path.join(ASSEMBLY_ROOTDIR, REFASM_SUBDIR)
REFASM_OUTPUT_DIR_PACKED = safe_format(os.path.join(ASSEMBLY_ROOTDIR, ASSEMBLY_SUBDIR), assembler = "refasm")
REFASM_LOG = os.path.join(REFASM_OUTPUT_DIR, "assembly.log")
REFASM_ASSEMBLED_CONTIGS = os.path.join(REFASM_OUTPUT_DIR, "contigs.fa")

# Evaluations

QUAST_ROOTDIR = os.path.join(EVALUATION_ROOTDIR, "quast-m{quast_mode}")
QUAST_SUBDIR = os.path.join(GENOME_ARGUMENT_STRING, GENOME_REFERENCE_ARGUMENT_STRING, GENOME_READS_ARGUMENT_STRING, ASSEMBLY_ARGUMENT_STRING)
QUAST_OUTPUT_DIR = os.path.join(QUAST_ROOTDIR, QUAST_SUBDIR)

RESOURCES_ROOTDIR = os.path.join(EVALUATION_ROOTDIR, "resources")
RESOURCES_SUBDIR = os.path.join(GENOME_ARGUMENT_STRING, GENOME_REFERENCE_ARGUMENT_STRING, GENOME_READS_ARGUMENT_STRING, ASSEMBLY_ARGUMENT_STRING)
RESOURCES_OUTPUT_DIR = os.path.join(RESOURCES_ROOTDIR, RESOURCES_SUBDIR)
RESOURCES_EVALUATION = os.path.join(RESOURCES_OUTPUT_DIR, "resources.json")

# Reports

REPORT_SUBDIR = os.path.join(REPORT_ROOTDIR, "{report_name}", "{report_file_name}")
REPORT_TEX = os.path.join(REPORT_SUBDIR, "report.tex")
REPORT_COMBINED_EAXMAX_PLOT = os.path.join(REPORT_SUBDIR, "combined_eaxmax_plot.pdf")
REPORT_NAME_FILE = os.path.join(REPORT_SUBDIR, "name.txt")
REPORT_HASHDIR = os.path.join(REPORT_ROOTDIR, "hashdir")
REPORT_PDF = os.path.join(REPORT_SUBDIR, "report.pdf")

AGGREGATED_REPORT_SUBDIR = os.path.join(REPORT_ROOTDIR, "{aggregated_report_name}")
AGGREGATED_REPORT_PDF = os.path.join(AGGREGATED_REPORT_SUBDIR, "aggregated_report.pdf")

UNIQUIFY_IDS_SCRIPT = "scripts/uniquify_fasta_ids.py"
CONVERT_VALIDATION_OUTPUTS_TO_LATEX_SCRIPT = "scripts/convert_validation_outputs_to_latex.py"
CREATE_AGGREGATED_WTDBG2_REPORT_SCRIPT = "scripts/create_aggregated_wtdbg2_report.py"
CREATE_COMBINED_EAXMAX_PLOT_SCRIPT = "scripts/create_combined_eaxmax_plot.py"
DOWNSAMPLE_FASTA_READS_SCRIPT = "scripts/downsample_fasta_reads.py"
FILTER_NW_FROM_REFERENCE_SCRIPT = "scripts/filter_nw_from_reference.py"
RETAIN_CM_IN_REFERENCE_SCRIPT = "scripts/retain_cm_in_reference.py"
FILTER_PLASMIDS_IN_REFERENCE_SCRIPT = "scripts/filter_plasmids_in_reference.py"
COMPUTE_GENOME_REFERENCE_LENGTH_SCRIPT = "scripts/compute_genome_reference_length.py"
WTDBG2_HODECO_SCRIPT = "scripts/wtdbg2_hodeco.py"

EXTERNAL_SOFTWARE_SCRIPTS_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "scripts")
RUST_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "rust_target")
IS_RUST_TESTED_MARKER = os.path.join(RUST_DIR, "is_rust_tested.marker")
IS_RUST_FETCHED_MARKER = os.path.join(RUST_DIR, "is_rust_fetched.marker")
RUST_BINARY = os.path.join(RUST_DIR, "release", "cli")
QUAST_BINARY = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "quast", "quast.py")
WTDBG2_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "wtdbg2")
WTDBG2_BINARY = os.path.join(WTDBG2_DIR, "wtdbg2")
WTDBG2_CONSENSUS_BINARY = os.path.join(WTDBG2_DIR, "wtpoa-cns")
FLYE_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "Flye")
FLYE_BINARY = os.path.join(FLYE_DIR, "bin", "flye")
SIM_IT_BINARY = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "sim-it", "sim-it.pl")
RATATOSK_BINARY = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "Ratatosk", "build", "src", "Ratatosk")
CONTIG_VALIDATOR_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "ContigValidator")
CONVERT_TO_GFA_BINARY = os.path.join(EXTERNAL_SOFTWARE_SCRIPTS_DIR, "convertToGFA.py")
SDSL_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "sdsl-lite")
MDBG_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "rust-mdbg")
MDBG_CARGO_TOML = os.path.join(MDBG_DIR, "Cargo.toml")
MDBG_BINARY = os.path.join(MDBG_DIR, "target", "release", "rust-mdbg")
MDBG_SIMPLIFY = os.path.join(MDBG_DIR, "utils", "magic_simplify")
MDBG_MULTI_K = os.path.join(MDBG_DIR, "utils", "multik")
LJA_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "LJA")
LJA_BINARY = os.path.join(LJA_DIR, "bin", "lja")
HIFIASM_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "hifiasm")
HIFIASM_BINARY = os.path.join(HIFIASM_DIR, "hifiasm")
HOMOPOLYMER_COMPRESS_RS_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "homopolymer-compress-rs")
HOMOPOLYMER_COMPRESS_RS_BINARY = os.path.join(HOMOPOLYMER_COMPRESS_RS_DIR, "target", "release", "homopolymer-compress")
WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "wtdbg2-homopolymer-decompression")
WTDBG2_HOMOPOLYMER_DECOMPRESSION_BINARY = os.path.join(WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR, "target", "release", "wtdbg2-homopolymer-decompression")
MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "minimap2-homopolymer-decompression")
MINIMAP2_HOMOPOLYMER_DECOMPRESSION_BINARY = os.path.join(MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR, "target", "release", "minimap2-homopolymer-decompression")
HISIM_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "HI.SIM")
HISIM_MODEL_BINARY = os.path.join(HISIM_DIR, "HImodel")
HISIM_SIM_BINARY = os.path.join(HISIM_DIR, "HIsim")
HISIM_FASTA_BINARY = os.path.join(HISIM_DIR, "HIfasta")
FASTK_DIR = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "FASTK")
FASTK_BINARY = os.path.join(FASTK_DIR, "FastK")
FASTK_SYMMEX_BINARY = os.path.join(FASTK_DIR, "Symmex")
FASTK_HISTEX_BINARY = os.path.join(FASTK_DIR, "Histex")

# TODO remove
ALGORITHM_PREFIX_FORMAT = os.path.join(DATADIR, "algorithms", "{arguments}")

################################
###### Cluster properties ######
################################

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
            if hasattr(wildcards, "read_source") and wildcards.read_source.startswith("hisim_"):
                return int(float(genome_properties["assembly_time_factor"]) * float(base_time_min)) * 2
            else:
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
        cluster = compute_genome_cluster_from_wildcards(wildcards, base_time_min, base_mem_mb)

        if cluster == "ukko":
            if mem >= 250_000:
                if time <= 1440 * 3:
                    return "bigmem,aurinko"
                else:
                    return "aurinko"
            elif time <= 60 * 8:
                return "short,medium,bigmem,aurinko"
            elif time <= 1440 * 2:
                return "medium,bigmem,aurinko"
            elif time <= 1440 * 3:
                return "bigmem,aurinko"
            else:
                return "aurinko" # aurinko can handle anything
            #else:
            #    sys.exit("No applicable queue for runtime " + str(time) + " (wildcards: " + str(wildcards) + ")")
        elif cluster == "kale":
            raise Exception("kale disabled")

            if mem >= 380_000:
                return "bigmem"
            elif time <= 1440:
                return "short,medium,long,bigmem"
            elif time <= 1440 * 7:
                return "medium,long,bigmem"
            elif time <= 1440 * 14:
                return "long,bigmem"
            else:
                sys.exit("No applicable queue for runtime " + str(time) + " (wildcards: " + str(wildcards) + ")")
        else:
            sys.exit("No applicable cluster for runtime " + str(time) + " (wildcards: " + str(wildcards) + ")")
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def compute_genome_cluster_from_wildcards(wildcards, base_time_min, base_mem_mb = 0):
    try:
        time = compute_genome_time_min_from_wildcards(wildcards, base_time_min)
        mem = compute_genome_mem_mb_from_wildcards(wildcards, base_mem_mb)

        #if time <= 1440 * 2:
        #    return "ukko"
        #elif time <= 1440 * 14:
        #    return "kale"
        #else:
        #    sys.exit("No applicable cluster for runtime " + str(time) + " (wildcards: " + str(wildcards) + ")")

        return "ukko" # run everything on ukko for now, using aurinko for long jobs
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

#################################
###### Global report rules ######
#################################

localrules: do_nothing
rule do_nothing:
    shell:  "echo 'No target specified'"

def get_all_report_files():
    try:
        result = []
        for report_name, report_definition in reports.items():
            for report_file_name, report_file_definition in report_definition["report_files"].items():
                report_file_arguments = report_file_definition.arguments
                #print(f"report_name: {report_name}")
                #print(f"report_file_arguments: {report_file_arguments}", flush = True)
                result.append(REPORT_PDF.format(report_name = report_name, report_file_name = str(report_file_arguments)))


        for aggregated_report_name in aggregated_reports.keys():
            result.append(AGGREGATED_REPORT_PDF.format(aggregated_report_name = aggregated_report_name))
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: report_all
rule report_all:
    input:  get_all_report_files(),
    threads: 1
    resources: mail_type = "END,FAIL,INVALID_DEPEND,REQUEUE"

def get_test_report_files():
    try:
        result = []
        for report_name, report_definition in reports.items():
            if not report_name.startswith("E.coli_HiFi"):
                continue

            for report_file_name, report_file_definition in report_definition["report_files"].items():
                report_file_arguments = report_file_definition.arguments
                #print(f"report_file_arguments: {report_file_arguments}", flush = True)
                result.append(REPORT_PDF.format(report_name = report_name, report_file_name = str(report_file_arguments)))


        for aggregated_report_name in aggregated_reports.keys():
            result.append(AGGREGATED_REPORT_PDF.format(aggregated_report_name = aggregated_report_name))
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: report_test
rule report_test:
    input:  get_test_report_files(),
    threads: 1
    resources: mail_type = "END,FAIL,INVALID_DEPEND,REQUEUE"

###############################
###### Report Generation ######
###############################

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
            #print(f"QUAST_OUTPUT_DIR: {QUAST_OUTPUT_DIR}")
            #print(f"column.arguments: {column.arguments}", flush = True)
            assembler = column.arguments.assembler_name()
            assembler_argument_string = ASSEMBLER_ARGUMENT_STRINGS[assembler]
            #print(assembler_argument_string, flush = True)

            assembler_arguments = column.arguments.assembler_arguments()
            #print(f"raw assembler_arguments: {assembler_arguments}", flush = True)
            assembler_arguments.setdefault("retain_cm", column.arguments["retain_cm"])
            #print(f"completed assembler_arguments: {assembler_arguments}", flush = True)
            assembler_arguments = assembler_argument_string.format(**assembler_arguments)
            #print(f"serialised assembler_arguments: {assembler_arguments}", flush = True)

            quasts.append(safe_format(QUAST_OUTPUT_DIR, assembler = assembler, assembler_arguments = assembler_arguments).format(**column.arguments))
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

def get_report_file_resources_evaluations_from_wildcards(wildcards):
    try:
        report_name = wildcards.report_name
        report_file_name = wildcards.report_file_name
        report_file = get_report_file(report_name, report_file_name)
        result = []
        for column in report_file.columns:
            #print(RESOURCES_EVALUATION)
            #print(column.arguments, flush = True)
            assembler = column.arguments.assembler_name()
            assembler_argument_string = ASSEMBLER_ARGUMENT_STRINGS[assembler]
            #print(assembler_argument_string, flush = True)
            #print(column.arguments.assembler_arguments(), flush = True)
            assembler_arguments = assembler_argument_string.format(**column.arguments.assembler_arguments())
            #print(assembler_arguments, flush = True)
            result.append(safe_format(RESOURCES_EVALUATION, assembler = assembler, assembler_arguments = assembler_arguments).format(**column.arguments))
        return result
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

            assembler = column.arguments.assembler_name()
            assembler_argument_string = ASSEMBLER_ARGUMENT_STRINGS[assembler]
            assembler_arguments = assembler_argument_string.format(**column.arguments.assembler_arguments())
            quast_output_dir = safe_format(QUAST_OUTPUT_DIR, assembler = assembler, assembler_arguments = assembler_arguments).format(**column.arguments)
            resources_evaluation_file = safe_format(RESOURCES_EVALUATION, assembler = assembler, assembler_arguments = assembler_arguments).format(**column.arguments)
            result += f"'{column.shortname}' '' '{quast_output_dir}' '' '{resources_evaluation_file}'"
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

rule create_single_report_tex:
    input:  quasts = get_report_file_quasts_from_wildcards,
            combined_eaxmax_plot = REPORT_COMBINED_EAXMAX_PLOT,
            resources_evaluations = get_report_file_resources_evaluations_from_wildcards,
            script = CONVERT_VALIDATION_OUTPUTS_TO_LATEX_SCRIPT,
    output: report = REPORT_TEX,
    params: genome_name = lambda wildcards: ", ".join(get_report_genome_names_from_wildcards(wildcards)),
            script_column_arguments = get_single_report_script_column_arguments_from_wildcards,
            name_file = REPORT_NAME_FILE,
            hashdir = REPORT_HASHDIR,
    wildcard_constraints:
            report_name = "[^/]+",
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        mkdir -p '{params.hashdir}'
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

rule create_aggregated_report_tex:
    input: source_reports = get_aggregated_report_file_source_report_paths_from_wildcards,
           script = "scripts/create_aggregated_wtdbg2_report.py",
    output: file = AGGREGATED_REPORT_PDF,
    params: source_reports_arg = lambda wildcards: "' '".join(get_aggregated_report_file_source_report_paths_from_wildcards(wildcards)),
            source_report_names_arg = lambda wildcards: "' '".join([report_name + "/" + report_file_name for report_name, report_file_name in iterate_aggregated_report_file_source_reports_short_names(wildcards.aggregated_report_name)]),            
    wildcard_constraints:
            report_name = "[^/]+",
    conda: "config/conda-latex-gen-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        python3 '{input.script}' --source-reports '{params.source_reports_arg}' --source-report-names '{params.source_report_names_arg}' --output '{output.file}'
        """

rule create_combined_eaxmax_graph:
    input:  quast_csvs = lambda wildcards: [os.path.join(q, "aligned_stats", "EAxmax_plot.csv") for q in get_report_file_quasts_from_wildcards(wildcards)],
            script = CREATE_COMBINED_EAXMAX_PLOT_SCRIPT,
    output: REPORT_COMBINED_EAXMAX_PLOT,
    params: input_quast_csvs = lambda wildcards, input: "' '".join([shortname + "' '" + quast for shortname, quast in zip(get_report_file_column_shortnames_from_wildcards(wildcards), input.quast_csvs)])
    wildcard_constraints:
            report_name = "[^/]+",
    conda:  "config/conda-seaborn-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        mkdir -p "$(dirname '{output}')"
        python3 '{input.script}' '{params.input_quast_csvs}' '{output}'
        """

rule png_to_pdf:
    input: "{file}.png"
    output: "{file}.image.pdf"
    conda: "config/conda-imagemagick-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: "convert {input} {output}"

rule latex:
    input: "{subpath}report.tex"
    output: "{subpath}report.pdf"
    conda: "config/conda-latex-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        tectonic '{input}'
        """

# def get_genome_reads_from_wildcards(wildcards):
#     try:
#         arguments = Arguments.from_str(wildcards.arguments)

#         if arguments.read_simulator_name() is None or arguments.read_simulator_name() == "none":
#             return GENOME_READS_FORMAT.format(genome = arguments["genome"])
#         else:
#             return GENOME_SIMULATED_READS_FASTA_FORMAT.format(genome = arguments["genome"], read_simulator_name = arguments.read_simulator_name(), read_simulator_arguments = arguments.read_simulator_arguments())
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# def get_genome_reference_from_wildcards(wildcards):
#     try:
#         arguments = Arguments.from_str(wildcards.arguments)
#         return GENOME_REFERENCE_FORMAT.format(genome = arguments["genome"])
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# rule find_wtdbg2_node_errors:
#     input:  nodes = os.path.join(WTDBG2_OUTPUT_DIR_PACKED, "wtdbg2.1.nodes"),
#             reference = GENOME_REFERENCE,
#             script = "scripts/find_wtdbg2_node_errors.py",
#     output: deviation_histogram = os.path.join(ASSEMBLY_OUTPUT_DIR, "wtdbg2_node_errors", "deviation_histogram.pdf"),
#     log:    log = os.path.join(ALGORITHM_PREFIX_FORMAT, "wtdbg2_node_errors", "wtdbg2_node_errors.log"),
#     params: output_prefix = os.path.join(ALGORITHM_PREFIX_FORMAT, "wtdbg2_node_errors") + "/",
#     conda:  "config/conda-wtdbg2-node-errors-env.yml"
#     threads: 1
#     shell:  "PYTHONUNBUFFERED=1 '{input.script}' '{input.nodes}' '{input.reference}' '{params.output_prefix}' 2>&1 | tee '{log.log}'"

########################
###### Injections ######
########################

def get_injectable_contigs_rust_cli_command_from_wildcards(wildcards):
    try:
        injection = wildcards.tig_injection

        if "unitigs" == injection:
            return "compute-unitigs"
        elif "trivial_omnitigs" == injection:
            return "compute-trivial-omnitigs --non-scc"
        elif "omnitigs" == injection:
            return "compute-omnitigs --linear-reduction --linear-reduction-use-scc"
        else:
            raise Exception(f"Wrong injection command in assembler arguments: {wildcards}")
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule compute_injectable_contigs_wtdbg2:
    input:  nodes = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.nodes"), tig_injection = "none"),
            reads = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.reads"), tig_injection = "none"),
            dot = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.dot"), tig_injection = "none"),
            raw_reads = GENOME_READS,
            binary = RUST_BINARY,
    output: file = os.path.join(WTDBG2_INJECTABLE_CONTIG_DIR, "contigwalks.ssv"),
            latex = os.path.join(WTDBG2_INJECTABLE_CONTIG_DIR, "compute_injectable_contigs.tex"),
    log:    log = os.path.join(WTDBG2_INJECTABLE_CONTIG_DIR, "compute_injectable_contigs.log"),
    params: command = get_injectable_contigs_rust_cli_command_from_wildcards,
    threads: 1
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 48_000),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 48_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 48_000),
    shell: "${{CONDA_PREFIX}}/bin/time -v '{input.binary}' --log-level Trace {params.command} --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{log.log}'"

def get_injectable_fragment_contigs_rust_cli_command_from_wildcards(wildcards):
    try:
        injection = wildcards.frg_injection

        if "unitigs" == injection:
            return "compute-unitigs"
        elif "trivial_omnitigs" == injection:
            return "compute-trivial-omnitigs --non-scc"
        elif "omnitigs" == injection:
            return "compute-omnitigs --linear-reduction --linear-reduction-use-scc"
        else:
            raise Exception("Wrong injection command in assembler arguments: {wildcards}")
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule compute_injectable_fragment_contigs_wtdbg2:
    input:  dot = lambda wildcards: safe_format(os.path.join(WTDBG2_OUTPUT_DIR, f"wtdbg2.{wildcards.frg_injection_stage}.dot"), frg_injection = "none", frg_injection_stage = "none"),
            binary = RUST_BINARY,
    output: file = os.path.join(WTDBG2_INJECTABLE_FRAGMENT_CONTIG_DIR, "contigwalks.ssv"),
            latex = os.path.join(WTDBG2_INJECTABLE_FRAGMENT_CONTIG_DIR, "compute_injectable_contigs.tex"),
    log:    log = os.path.join(WTDBG2_INJECTABLE_FRAGMENT_CONTIG_DIR, "compute_injectable_contigs.log"),
    params: command = get_injectable_fragment_contigs_rust_cli_command_from_wildcards,
    threads: 1
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 48_000),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 48_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 48_000),
    shell: "${{CONDA_PREFIX}}/bin/time -v '{input.binary}' --log-level Trace {params.command} --file-format dot --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{log.log}'"

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

def get_raw_gfa_assembly_file_from_wildcards(wildcards):
    try:
        arguments = Arguments.from_str(wildcards.arguments)
        assembler_name = arguments.assembler_name()

        if assembler_name == "hifiasm":
            result = os.path.join(HIFIASM_PREFIX_FORMAT, "hifiasm", "assembly.r_utg.gfa")
        else:
            raise Exception("Assembler {} does not support gfa (in arguments: {})".format(assembler_name, arguments))

        arguments.retain_raw_assembly_arguments()
        return result.format(genome = arguments.genome(), arguments = str(arguments))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

localrules: select_gfa_assembler
rule select_gfa_assembler:
    input:  raw_assembly_from_assembler = get_raw_gfa_assembly_file_from_wildcards,
    output: raw_assembly = ALGORITHM_PREFIX_FORMAT + "raw_assembly.gfa",
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
        elif postprocessor_name == "gfa_trivial_omnitigs":
            return os.path.join(ALGORITHM_PREFIX_FORMAT, "gfa_trivial_omnitigs/trivial_omnitigs.fa")
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

rule run_contigbreaker:
    input:  contigs = ALGORITHM_PREFIX_FORMAT + "raw_assembly.fa",
            reads = GENOME_READS,
            script = "tools/contigbreaker/contigbreaker.py"
    output: broken_contigs = ALGORITHM_PREFIX_FORMAT + "contigbreaker/broken_contigs.fa",
            completed = touch(ALGORITHM_PREFIX_FORMAT + "contigbreaker/broken_contigs.fa.completed"),
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 60_000),
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               cpus = MAX_THREADS,
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 60_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 60_000),
    conda: "tools/contigbreaker/environment.yml"
    shell: "'{input.script}' --threads {threads} --input-contigs '{input.contigs}' --input-reads '{input.reads}' --output-contigs '{output.broken_contigs}'"

rule run_gfa_trivial_omnitigs:
    input:  contigs = os.path.join(ALGORITHM_PREFIX_FORMAT, "raw_assembly.gfa"),
            binary = RUST_BINARY,
    output: trivial_omnitigs = os.path.join(ALGORITHM_PREFIX_FORMAT, "gfa_trivial_omnitigs", "trivial_omnitigs.fa"),
    log:    log = os.path.join(ALGORITHM_PREFIX_FORMAT, "gfa_trivial_omnitigs", "trivial_omnitigs.log"),
    threads: 1
    resources:
            mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 10_000),
            time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
            cpus = 1,
            queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 10_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 10_000),
    shell: "'{input.binary}' compute-trivial-omnitigs --non-scc --file-format hifiasm --input '{input.contigs}' --output '{output.trivial_omnitigs}' 2>&1 | tee '{log.log}'"


####################
###### wtdbg2 ######
####################

rule wtdbg2:
    input:  reads = GENOME_READS,
            reference_length = UNFILTERED_GENOME_REFERENCE_LENGTH,
            binary = WTDBG2_BINARY,
    output: original_nodes = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.1.nodes"),
            nodes = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.nodes"),
            reads = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.reads"),
            dot = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.dot.gz"),
            ctg_dot = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.dot.gz"),
            clips = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.clps"),
            kbm = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.kbm"),
            ctg_lay = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay.gz"),
            frg_dot = [os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.{}.frg.dot.gz".format(i)) for i in range(1, 11)],
    log:    log = WTDBG2_LOG,
    params: output_prefix = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2"),
            fragment_correction_steps = lambda wildcards: f"--fragment-correction-steps {wildcards.fragment_correction_steps}" if wildcards.fragment_correction_steps != "all" else "",
    wildcard_constraints:
            tig_injection = "none",
            frg_injection = "none",
            frg_injection_stage = "none",
            skip_fragment_assembly = "no",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 100_000),
    shell: """
        read -r REFERENCE_LENGTH < '{input.reference_length}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -x {wildcards.wtdbg2_mode} -g $REFERENCE_LENGTH -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --dump-kbm '{output.kbm}' {params.fragment_correction_steps} 2>&1 | tee '{log.log}'
    """

rule wtdbg2_skip_fragment_assembly:
    input:  reads = GENOME_READS,
            reference_length = UNFILTERED_GENOME_REFERENCE_LENGTH,
            binary = WTDBG2_BINARY,
    output: original_nodes = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.1.nodes"),
            nodes = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.nodes"),
            reads = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.reads"),
            dot = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.3.dot.gz"),
            ctg_dot = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.dot.gz"),
            clips = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.clps"),
            kbm = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.kbm"),
            ctg_lay = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay.gz"),
    log:    log = WTDBG2_LOG,
    params: output_prefix = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2"),
            fragment_correction_steps = lambda wildcards: f"--fragment-correction-steps {wildcards.fragment_correction_steps}" if wildcards.fragment_correction_steps != "all" else "",
    wildcard_constraints:
            tig_injection = "none",
            frg_injection = "none",
            frg_injection_stage = "none",
            skip_fragment_assembly = "yes",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 100_000),
    shell: """
        read -r REFERENCE_LENGTH < '{input.reference_length}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -x {wildcards.wtdbg2_mode} -g $REFERENCE_LENGTH -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --dump-kbm '{output.kbm}' --skip-fragment-assembly {params.fragment_correction_steps} 2>&1 | tee '{log.log}'
    """

rule wtdbg2_with_injected_contigs:
    input:  reads = GENOME_READS,
            contigs = os.path.join(WTDBG2_INJECTABLE_CONTIG_DIR, "contigwalks.ssv"),
            #fragment_contigs = get_wtdbg2_injectable_fragment_contigs_from_wildcards,
            cached_nodes = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.1.nodes"), tig_injection = "none", fragment_correction_steps = "all"),
            cached_clips = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.clps"), tig_injection = "none", fragment_correction_steps = "all"),
            cached_kbm = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.kbm"), tig_injection = "none", fragment_correction_steps = "all"),
            reference_length = UNFILTERED_GENOME_REFERENCE_LENGTH,
            edge_cov = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.edge-cov"), tig_injection = "none", fragment_correction_steps = "all"),
            binary = WTDBG2_BINARY,
    output: ctg_lay = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay.gz"),
    log:    log = WTDBG2_LOG,
    params: output_prefix = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2"),
            skip_fragment_assembly = lambda wildcards: "--skip-fragment-assembly" if wildcards.skip_fragment_assembly == "yes" else "",
            fragment_correction_steps = lambda wildcards: f"--fragment-correction-steps {wildcards.fragment_correction_steps}" if wildcards.fragment_correction_steps != "all" else "",
    wildcard_constraints:
            tig_injection = "((?!none).)*",
            frg_injection = "none",
            frg_injection_stage = "none",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 100_000),
    shell: """
        read -r REFERENCE_LENGTH < '{input.reference_length}'
        read -r EDGE_COV < '{input.edge_cov}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -x {wildcards.wtdbg2_mode} -g $REFERENCE_LENGTH -e $EDGE_COV -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --load-nodes '{input.cached_nodes}' --load-clips '{input.cached_clips}' --load-kbm '{input.cached_kbm}' --inject-unitigs '{input.contigs}' {params.skip_fragment_assembly} {params.fragment_correction_steps} 2>&1 | tee '{log.log}'
    """

rule wtdbg2_with_injected_fragment_contigs:
    input:  reads = GENOME_READS,
            fragment_contigs = os.path.join(WTDBG2_INJECTABLE_FRAGMENT_CONTIG_DIR, "contigwalks.ssv"),
            cached_nodes = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.1.nodes"), frg_injection = "none", frg_injection_stage = "none"),
            cached_clips = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.clps"), frg_injection = "none", frg_injection_stage = "none"),
            cached_kbm = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.kbm"), frg_injection = "none", frg_injection_stage = "none"),
            reference_length = UNFILTERED_GENOME_REFERENCE_LENGTH,
            edge_cov = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.edge-cov"), frg_injection = "none", frg_injection_stage = "none"),
            binary = WTDBG2_BINARY,
    output: ctg_lay = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay.gz"),
    log:    log = WTDBG2_LOG,
    params: output_prefix = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2"),
    wildcard_constraints:
            tig_injection = "none",
            frg_injection = "((?!none).)*",
            frg_injection_stage = "((?!none).)*",
            skip_fragment_assembly = "no",
            fragment_correction_steps = "all",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 100_000),
    shell: """
        read -r REFERENCE_LENGTH < '{input.reference_length}'
        read -r EDGE_COV < '{input.edge_cov}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -x {wildcards.wtdbg2_mode} -g $REFERENCE_LENGTH -e $EDGE_COV -i '{input.reads}' -t {threads} -fo '{params.output_prefix}' --load-nodes '{input.cached_nodes}' --load-clips '{input.cached_clips}' --load-kbm '{input.cached_kbm}' --inject-fragment-unitigs '{input.fragment_contigs}' 2>&1 | tee '{log.log}'
    """

rule wtdbg2_extract:
    input:  file = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.{subfile}.gz"),
    output: file = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.{subfile}"),
    log:    WTDBG2_EXTRACT_LOG,
    params: working_directory = lambda wildcards, input: os.path.dirname(input.file),
            abslog = lambda wildcards: os.path.abspath(WTDBG2_EXTRACT_LOG.format(**wildcards)),
    wildcard_constraints:
            hodeco_consensus = "none",
    conda:  "config/conda-extract-env.yml"
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 360, 1_000),
    shell:  """
        cd '{params.working_directory}'
        ${{CONDA_PREFIX}}/bin/time -v gunzip -k wtdbg2.{wildcards.subfile}.gz 2>&1 | tee '{params.abslog}'
        """

rule wtdbg2_consensus:
    input:  contigs = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay"),
            binary = WTDBG2_CONSENSUS_BINARY,
    output: consensus = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.raw.fa"),
    log:    log = WTDBG2_CONSENSUS_LOG,
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 8_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360, 8_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 360, 8_000),
    shell: "${{CONDA_PREFIX}}/bin/time -v {input.binary} -t {threads} -i '{input.contigs}' -fo '{output.consensus}' 2>&1 | tee '{log.log}'"

rule wtdbg2_transform_ctg_lay_hodeco_simple:
    input:  contigs = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay"), hodeco_consensus = "none", homopolymer_compression = "yes"),
            normal_reads = safe_format(GENOME_READS, homopolymer_compression = "none"),
            #hoco_reads = safe_format(GENOME_READS, homopolymer_compression = "yes"),
            binary = WTDBG2_HOMOPOLYMER_DECOMPRESSION_BINARY,
    output: contigs = os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.ctg.lay"),
    log:    log = os.path.join(WTDBG2_OUTPUT_DIR, "transform_ctg_lay_hodeco_simple.log"),
    #conda:  "config/conda-biopython-env.yml"
    wildcard_constraints:
            hodeco_consensus = "simple",
            homopolymer_compression = "none",
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 18_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360, 18_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 360, 18_000),
    shell:  "${{CONDA_PREFIX}}/bin/time -v {input.binary} --input {input.contigs} --output {output.contigs} --normal-reads {input.normal_reads} --compute-threads {threads} 2>&1 | tee '{log.log}'"


localrules: link_wtdbg2_contigs
rule link_wtdbg2_contigs:
    input:  contigs = os.path.join(WTDBG2_OUTPUT_DIR_PACKED, "wtdbg2.raw.fa"),
    output: contigs = ASSEMBLED_CONTIGS,
    wildcard_constraints:
            assembler = "wtdbg2",
    shell:  "ln -sr -T '{input.contigs}' '{output.contigs}'"

rule wtdbg2_find_edge_cov:
    input:  WTDBG2_LOG,
    output: os.path.join(WTDBG2_OUTPUT_DIR, "wtdbg2.edge-cov"),
    shell:  """
        grep 'Set --edge-cov to ' '{input}' | sed 's/.*Set --edge-cov to //g' > '{output}'
        """

##################
###### Flye ######
##################

rule flye:
    input:  reads = GENOME_READS,
            reference_length = UNFILTERED_GENOME_REFERENCE_LENGTH,
            script = FLYE_BINARY,
    output: contigs = os.path.join(FLYE_OUTPUT_DIR, "flye", "assembly.fasta"),
            directory = directory(os.path.join(FLYE_OUTPUT_DIR, "flye")),
    log:    log = FLYE_LOG,
    params: output_directory = os.path.join(FLYE_OUTPUT_DIR, "flye"),
    #conda: "config/conda-flye-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 250_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440, 250_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 1440, 250_000),
    shell:  """
        read -r REFERENCE_LENGTH < '{input.reference_length}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.script}' -g $REFERENCE_LENGTH -t {threads} -o '{params.output_directory}' --{wildcards.flye_mode} '{input.reads}' 2>&1 | tee '{log.log}'
        """

localrules: link_flye_contigs
rule link_flye_contigs:
    input:  contigs = os.path.join(FLYE_OUTPUT_DIR_PACKED, "flye", "assembly.fasta"),
    output: contigs = ASSEMBLED_CONTIGS,
    wildcard_constraints:
            assembler = "flye",
    shell:  "ln -sr -T '{input.contigs}' '{output.contigs}'"

#####################
###### Hifiasm ######
#####################

rule hifiasm:
    input:  reads = GENOME_READS,
            binary = HIFIASM_BINARY,
    output: contigs = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.p_ctg.gfa"),
            alternate_contigs = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.a_ctg.gfa"),
            raw_unitigs = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.r_utg.gfa"),
            raw_unitigs_bubpop = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.bubpop.r_utg.gfa"),
            primary_unitigs = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.p_utg.gfa"),
            directory = directory(os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm")),
    params: output_prefix = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly"),
    log:    HIFIASM_LOG
    wildcard_constraints:
            contig_algorithm = "builtin",
    conda:  "config/conda-hifiasm-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 50_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 50_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 50_000),
    shell: "${{CONDA_PREFIX}}/bin/time -v '{input.binary}' --primary -t {threads} -o '{params.output_prefix}' '{input.reads}' 2>&1 | tee '{log}'"

rule hifiasm_gfa_to_fa:
    input:  gfa = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.p_ctg.gfa"),
            alternate_gfa = os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.a_ctg.gfa"),
    output: fa = HIFIASM_ASSEMBLED_CONTIGS,
    wildcard_constraints:
            contig_algorithm = "builtin",
    wildcard_constraints:
            assembler = "hifiasm",
    run:
            with open(input.gfa, 'r') as input_file, open(input.alternate_gfa, 'r') as alternate_input_file, open(output.fa, 'w') as output_file:
                for line in itertools.chain(input_file, alternate_input_file):
                    if line[0] != "S":
                        continue

                    columns = line.split("\t")
                    print(f"Writing contig {columns[1]}...")
                    output_file.write(">{}\n{}\n".format(columns[1], columns[2]))
                print(f"Wrote all contigs")

rule hifiasm_trivial_omnitigs:
    input:  unitigs = lambda wildcards: safe_format(os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.{graph_variant}.gfa"), contig_algorithm = "builtin", graph_variant = wildcards.contig_algorithm[len("trivial_omnitigs_"):]),
            rust_binary = RUST_BINARY,
    output: contigs = HIFIASM_ASSEMBLED_CONTIGS,
            latex = os.path.join(HIFIASM_OUTPUT_DIR, "compute_injectable_contigs.tex"),
    wildcard_constraints:
            contig_algorithm = r"trivial_omnitigs_[a-z_\.]+",
    log:    log = os.path.join(HIFIASM_OUTPUT_DIR, "compute_injectable_contigs.log"),
    shell: "${{CONDA_PREFIX}}/bin/time -v '{input.rust_binary}' compute-trivial-omnitigs --file-format hifiasm --input '{input.unitigs}' --output '{output.contigs}' --latex '{output.latex}' --non-scc 2>&1 | tee '{log.log}'"

rule hifiasm_omnitigs:
    input:  unitigs = lambda wildcards: safe_format(os.path.join(HIFIASM_OUTPUT_DIR, "hifiasm", "assembly.{graph_variant}.gfa"), contig_algorithm = "builtin", graph_variant = wildcards.contig_algorithm[len("omnitigs_"):]),
            rust_binary = RUST_BINARY,
    output: contigs = HIFIASM_ASSEMBLED_CONTIGS,
            latex = os.path.join(HIFIASM_OUTPUT_DIR, "compute_injectable_contigs.tex"),
    wildcard_constraints:
            contig_algorithm = r"omnitigs_[a-z_\.]+",
    log:    log = os.path.join(HIFIASM_OUTPUT_DIR, "compute_injectable_contigs.log"),
    shell: "${{CONDA_PREFIX}}/bin/time -v '{input.rust_binary}' compute-omnitigs --file-format hifiasm --input '{input.unitigs}' --output '{output.contigs}' --latex '{output.latex}' --linear-reduction --linear-reduction-use-scc 2>&1 | tee '{log.log}'"

##################
###### mdbg ######
##################

rule mdbg_multik:
    input:  reads = GENOME_SINGLE_LINE_READS,
            script = MDBG_MULTI_K,
            binary = MDBG_BINARY,
    output: contigs = MDBG_ASSEMBLED_CONTIGS,
    params: output_prefix = os.path.join(MDBG_OUTPUT_DIR, "contigs"),
            original_contigs = os.path.join(MDBG_OUTPUT_DIR, "contigs-final.msimpl.fa"),
    wildcard_constraints:
            mdbg_mode = "multik",
    log:    log = MDBG_LOG,
    conda:  "config/conda-mdbg-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 1440, 100_000),
    shell:  """
        RUST_BACKTRACE=full ${{CONDA_PREFIX}}/bin/time -v '{input.script}' '{input.reads}' '{params.output_prefix}' {threads} 2>&1 | tee '{log.log}'
        ln -sr -T '{params.original_contigs}' '{output.contigs}'
        """

rule mdbg_D_melanogaster:
    input:  reads = GENOME_SINGLE_LINE_READS,
            simplify_script = MDBG_SIMPLIFY,
            binary = MDBG_BINARY,
    output: contigs = MDBG_ASSEMBLED_CONTIGS,
    params: output_prefix = os.path.join(MDBG_OUTPUT_DIR, "contigs"),
            original_contigs = os.path.join(MDBG_OUTPUT_DIR, "contigs.msimpl.fa"),
    wildcard_constraints:
            mdbg_mode = "D_melanogaster",
    log:    log = MDBG_LOG,
    conda:  "config/conda-mdbg-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 1440, 100_000),
    shell:  """
        RUST_BACKTRACE=full ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' '{input.reads}' -k 35 -l 12 --density 0.002 --threads {threads} --prefix '{params.output_prefix}' 2>&1 | tee '{log.log}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.simplify_script}' '{params.output_prefix}' 2>&1 | tee -a '{log.log}'
        ln -sr -T '{params.original_contigs}' '{output.contigs}'
        """

rule mdbg_HG002:
    input:  reads = GENOME_SINGLE_LINE_READS,
            simplify_script = MDBG_SIMPLIFY,
            binary = MDBG_BINARY,
    output: contigs = MDBG_ASSEMBLED_CONTIGS,
    params: output_prefix = os.path.join(MDBG_OUTPUT_DIR, "contigs"),
            original_contigs = os.path.join(MDBG_OUTPUT_DIR, "contigs.msimpl.fa"),
    wildcard_constraints:
            mdbg_mode = "HG002",
    log:    log = MDBG_LOG,
    conda:  "config/conda-mdbg-env.yml"
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 1440, 100_000),
    shell:  """
        RUST_BACKTRACE=full ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' '{input.reads}' -k 21 -l 14 --density 0.003 --threads {threads} --prefix '{params.output_prefix}' 2>&1 | tee '{log.log}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.simplify_script}' '{params.output_prefix}' 2>&1 | tee -a '{log.log}'
        ln -sr -T '{params.original_contigs}' '{output.contigs}'
        """

#################
###### LJA ######
#################

rule lja:
    input:  reads = GENOME_READS,
            binary = LJA_BINARY,
    output: contigs = LJA_ASSEMBLED_CONTIGS,
    params: output_dir = os.path.join(LJA_OUTPUT_DIR, "output"),
            original_contigs = os.path.join(LJA_OUTPUT_DIR, "output", "assembly.fasta"),
    log:    log = LJA_LOG,
    threads: MAX_THREADS
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 1440),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 1440, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 1440, 100_000),
    shell:  """
        mkdir -p '{params.output_dir}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -t {threads} -o '{params.output_dir}' --reads '{input.reads}' 2>&1 | tee '{log.log}'
        ln -sr -T '{params.original_contigs}' '{output.contigs}'
        """

####################
###### HiCanu ######
####################

rule hicanu:
    input:  reads = GENOME_READS,
            reference_length = UNFILTERED_GENOME_REFERENCE_LENGTH,
    output: contigs = CANU_ASSEMBLED_CONTIGS,
    params: output_dir = os.path.join(CANU_OUTPUT_DIR, "output"),
            original_contigs = os.path.join(CANU_OUTPUT_DIR, "output", "assembly.contigs.fasta"),
    log:    log = CANU_LOG,
    threads: MAX_THREADS,
    conda:  "config/conda-canu-env.yml"
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 720),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 720, 100_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 720, 100_000),
    shell:  """
        read -r REFERENCE_LENGTH < '{input.reference_length}'
        mkdir -p '{params.output_dir}'
        ${{CONDA_PREFIX}}/bin/time -v canu -assemble -p assembly -d '{params.output_dir}' genomeSize=$REFERENCE_LENGTH useGrid=false -pacbio-hifi '{input.reads}' 2>&1 | tee '{log.log}'
        ln -sr -T '{params.original_contigs}' '{output.contigs}'
        """

####################
###### RefAsm ######
####################

localrules: refasm
rule refasm:
    input:  reference = UNFILTERED_GENOME_REFERENCE,
    output: contigs = REFASM_ASSEMBLED_CONTIGS,
    log:    log = REFASM_LOG,
    threads: 1,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 10),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 10, 100),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 10, 100),
    shell:  """
        ln -sr -T '{input.reference}' '{output.contigs}'
        """

#############################
###### Read Simulation ######
#############################

# def get_read_simulator_args_from_wildcards(wildcards):
#     try:
#         read_simulator_arguments = Arguments.from_str(wildcards.read_simulator_arguments)
#         cli_arguments = read_simulator_arguments.get("cli_arguments", None)
#         if cli_arguments is None:
#             return ""

#         return cli_arguments.to_argument_string()
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# rule simulate_perfect_reads:
#     input:  reference = GENOME_REFERENCE,
#             script = "scripts/simulate_perfect_reads.py",
#     output: simulated_reads = safe_format(GENOME_SIMULATED_READS_FASTA_FORMAT, read_simulator_name = "perfect"),
#     log:    log = safe_format(GENOME_SIMULATED_READS_FASTA_FORMAT, read_simulator_name = "perfect") + ".log",
#     params: cli_arguments = get_read_simulator_args_from_wildcards,
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
#         mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1000),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 1000),
#     conda:  "config/conda-simulate-perfect-reads-env.yml",
#     shell:  "'{input.script}' {params.cli_arguments} --reference '{input.reference}' --output '{output.simulated_reads}'"

# rule simulate_hifi_reads_bbmap:
#     input:  reference = GENOME_REFERENCE,
#     output: simulated_reads = safe_format(GENOME_SIMULATED_READS_FASTA_FORMAT, read_simulator_name = "bbmap_hifi"),
#     log:    log = safe_format(GENOME_SIMULATED_READS_FASTA_FORMAT, read_simulator_name = "bbmap_hifi") + ".log",
#     params: working_directory = lambda wildcards, output: os.path.dirname(output.simulated_reads),
#             mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
#         mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 4000),
#     conda:  "config/conda-bbmap-env.yml"
#     shell:  """
#         REFERENCE=$(realpath -s '{input.reference}')
#         OUTPUT=$(realpath -s '{output.simulated_reads}')
#         LOG=$(realpath -s '{log.log}')

#         cd '{params.working_directory}'
#         randomreads.sh build=1 \
#         -Xmx{params.mem_mb}m \
#         ow=t seed=1 \
#         ref="$REFERENCE" \
#         simplenames=t \
#         pacbio=t pbmin=0.001 pbmax=0.01 \
#         coverage=30 paired=f \
#         gaussianlength=t \
#         minlength=9000 midlength=10000 maxlength=12000 \
#         out="$OUTPUT" 2>&1 | tee "$LOG"
#         """

# rule simulate_ccs_reads_simlord:
#     input:  reference = GENOME_REFERENCE,
#     output: simulated_reads = safe_format(GENOME_SIMULATED_READS_FASTQ_FORMAT, read_simulator_name = "simlord"),
#     log:    log = safe_format(GENOME_SIMULATED_READS_FASTQ_FORMAT, read_simulator_name = "simlord") + ".log",
#     params: working_directory = lambda wildcards, output: os.path.dirname(output.simulated_reads),
#             mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
#             cli_arguments = get_read_simulator_args_from_wildcards,
#             output_prefix = safe_format(GENOME_SIMULATED_READS_PREFIX_FORMAT, read_simulator_name = "simlord"),
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
#         mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 4000),
#     conda:  "config/conda-simlord-env.yml"
#     shell:  """
#         simlord --read-reference '{input.reference}' --no-sam {params.cli_arguments} '{params.output_prefix}'
#         mv '{params.output_prefix}.fastq' '{output.simulated_reads}'
#     """

# rule fastq_to_fasta:
#     input:  fastq = os.path.join(DATADIR, "{path}.fq"),
#     output: fasta = os.path.join(DATADIR, "{path}.fa"),
#     shell:  "awk '(NR-1)%4<2' '{input.fastq}' > '{output.fasta}'"

# rule simulate_hifi_reads_sim_it:
#     input:  reference = GENOME_REFERENCE,
#             script = SIM_IT_BINARY,
#     output: simulated_reads = safe_format(GENOME_SIMULATED_READS_FASTA_FORMAT, read_simulator_name = "sim-it_hifi"),
#     log:    log = safe_format(GENOME_SIMULATED_READS_FASTA_FORMAT, read_simulator_name = "sim-it_hifi") + ".log",
#     params: working_directory = lambda wildcards, output: os.path.dirname(output.simulated_reads),
#             mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
#             cli_arguments = get_read_simulator_args_from_wildcards,
#             output_prefix = safe_format(GENOME_SIMULATED_READS_PREFIX_FORMAT, read_simulator_name = "sim-it_hifi"),
#             config_file = safe_format(GENOME_SIMULATED_READS_PREFIX_FORMAT, read_simulator_name = "sim-it_hifi") + ".config",
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
#         mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 4000),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 4000),
#     conda:  "config/conda-perl-env.yml"
#     shell:  """
#         echo "Project:
# -----------------------------
# Project name             = Test
# Reference sequence       = {input.reference}
# Replace ambiguous nts(N) = 


# Structural variation:
# -----------------------------
# VCF input                = 
# Foreign sequences        =

# Deletions                = 0
# Length (bp)              = 30-150000

# Insertions               = 0
# Length (bp)              = 30-100000

# Tandem duplications      = 0
# Length (bp)              = 50-10000
# Copies                   = 1-20

# Inversions               = 0
# Length (bp)              = 150-1300000

# Complex substitutions    = 0
# Length (bp)              = 30-100000

# Inverted duplications    = 0
# Length (bp)              = 150-350000

# Heterozygosity           = 60%


# Long Read simulation:
# -----------------------------
# Coverage                 = 15
# Median length            = 10000
# Length range             = 7000-13000
# Accuracy                 = 99%
# Error profile            = external-software/sim-it/error_profile_PB_Sequel_CCS_hifi.txt" > '{params.config_file}'

#         mkdir -p '{params.output_prefix}'
#         perl '{input.script}' -c '{params.config_file}' -o '{params.output_prefix}' 2>&1 | tee '{log.log}'
#         ln -sr -T '{params.output_prefix}/Nanopore_Test.fasta' '{output.simulated_reads}'
#     """


#########################################
###### Long Read Input Preparation ######
#########################################

# def url_file_format(url):
#     parts = url.split('.')
#     if parts[-1] == "gz":
#         return parts[-2]
#     else:
#         return parts[-1]

# def read_url_file_format(genome):
#     try:
#         format = url_file_format(genomes[genome]["urls"][0])

#         if format == "fasta":
#             return "fa"
#         elif genomes[genome].get("format", None) == "sra":
#             return "sra"
#         else:
#             return format
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# ruleorder: download_raw_source_reads > download_packed_source_reads > convert_source_reads > extract

# localrules: download_raw_source_reads
# rule download_raw_source_reads:
#     output: file = DATADIR + "downloads/{genome}/reads-{index}/reads.{format}",
#             completed = touch(DATADIR + "downloads/{genome}/reads-{index}/reads.{format}.completed"),
#     params: url = lambda wildcards: genomes[wildcards.genome]["urls"][int(wildcards.index)],
#             url_format = lambda wildcards: read_url_file_format(wildcards.genome),
#     wildcard_constraints:
#         format = "(fa|bam|sra)",
#         index = "\d+",
#         genome = "((?!downloads).)*",
#     conda: "config/conda-download-env.yml"
#     threads: 1
#     shell: """
#         mkdir -p "$(dirname '{output.file}')"

#         if [ '{params.url_format}' != '{wildcards.format}' ]; then
#             echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
#             exit 1
#         fi

#         wget --progress=dot:mega -O '{output.file}' '{params.url}'
#         """

# localrules: download_packed_source_reads
# rule download_packed_source_reads:
#     output: file = DATADIR + "downloads/{genome}/packed-reads-{index}/reads.{format}.gz",
#             completed = touch(DATADIR + "downloads/{genome}/packed-reads-{index}/reads.{format}.gz.completed"),
#     params: url = lambda wildcards: genomes[wildcards.genome]["urls"][int(wildcards.index)],
#             url_format = lambda wildcards: read_url_file_format(wildcards.genome),
#             checksum = lambda wildcards: genomes[wildcards.genome]["checksums"][int(wildcards.index)] if "checksums" in genomes[wildcards.genome] else "",
#     wildcard_constraints:
#         format = "(fa|bam|sra)",
#         index = "\d+",
#         genome = "((?!downloads).)*",
#     conda: "config/conda-download-env.yml"
#     threads: 1
#     shell: """
#         mkdir -p "$(dirname '{output.file}')"

#         if [ '{params.url_format}' != '{wildcards.format}' ]; then
#             echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
#             exit 1
#         fi

#         wget --progress=dot:mega -O '{output.file}' '{params.url}'

#         if [ ! -z "{params.checksum}" ]; then
#             md5sum -c <<< "{params.checksum} {output.file}"
#         else
#             echo "Assuming file '{output.file}' was downloaded correctly since no checksum was provided."
#         fi
#         """

# def read_raw_input_file_name(wildcards):
#     try:
#         genome_name = wildcards.genome
#         genome_properties = genomes[genome_name]

#         if genome_properties["urls"][0].split('.')[-1] == "gz":
#             input_file_name = DATADIR + "downloads/" + genome_name + "/packed-reads-" + wildcards.index + "/reads." + read_url_file_format(wildcards.genome)
#         else:
#             input_file_name = DATADIR + "downloads/" + genome_name + "/reads-" + wildcards.index + "/reads." + read_url_file_format(wildcards.genome)
#         return input_file_name
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# rule convert_source_reads:
#     input: file = read_raw_input_file_name,
#     output: file = DATADIR + "downloads/{genome}/reads-{index}/reads.converted.fa",
#             completed = touch(DATADIR + "downloads/{genome}/reads-{index}/reads.converted.fa.completed"),
#     params: file_format = lambda wildcards: read_url_file_format(wildcards.genome)
#     wildcard_constraints:
#         index = "\d+",
#         genome = "((?!downloads).)*",
#     conda: "config/conda-convert-reads-env.yml"
#     threads: 1
#     shell: """
#         if [ '{params.file_format}' == 'bam' ]; then
#             samtools fasta '{input.file}' > '{output.file}'
#         elif [ '{params.file_format}' == 'sra' ]; then
#             fastq-dump --stdout --fasta default '{input.file}' > '{output.file}'
#         else
#             ln -sr -T '{input.file}' '{output.file}'
#         fi
#         """

# rule combine_reads:
#     input: files = lambda wildcards: expand(DATADIR + "downloads/{{genome}}/reads-{index}/reads.converted.fa", index=range(len(genomes[wildcards.genome]["urls"])))
#     output: reads = DATADIR + "downloads/{genome}/reads/raw_reads.fa",
#             completed = touch(DATADIR + "downloads/{genome}/reads/raw_reads.fa.completed"),
#     params: input_list = lambda wildcards, input: "'" + "' '".join(input.files) + "'"
#     wildcard_constraints:
#         genome = "((?!downloads).)*",
#     threads: 1
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60),
#     shell: "cat {params.input_list} > '{output.reads}'"

# def read_input_file_name(wildcards):
#     try:
#         genome_name = wildcards.genome
#         if genome_name in genomes:
#             return DATADIR + "downloads/" + genome_name + "/reads/raw_reads.uniquified.fa"
#         elif genome_name in corrected_genomes:
#             return DATADIR + "corrected_reads/" + genome_name + "/corrected_reads.fa"
#         else:
#             sys.exit("genome name not found: " + genome_name)
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# localrules: link_reads
# rule link_reads:
#     input: file = read_input_file_name,
#     #input: file = DATADIR + "corrected_reads/{genome}/corrected_reads.fa"
#     output: file = GENOME_READS_FORMAT,
#     wildcard_constraints:
#         genome = "((?!downloads).)*",
#     threads: 1
#     shell: "ln -sr -T '{input.file}' '{output.file}'"

# localrules: download_reference_raw
# rule download_reference_raw:
#     output: reference = DATADIR + "downloads/{genome}/reference/raw_reference.fa",
#             completed = touch(DATADIR + "downloads/{genome}/reference/raw_reference.fa.completed"),
#     params: url = lambda wildcards, output: genomes[wildcards.genome]["reference"]
#     wildcard_constraints:
#         genome = "((?!downloads).)*",
#     conda: "config/conda-download-env.yml"
#     threads: 1
#     shell: """
#         mkdir -p "$(dirname '{output.reference}')"
#         wget --progress=dot:mega -O '{output.reference}' '{params.url}'
#         """

# localrules: download_reference_gzip
# rule download_reference_gzip:
#     output: reference = DATADIR + "downloads/{genome}/packed-reference/raw_reference.fa.gz",
#             completed = touch(DATADIR + "downloads/{genome}/packed-reference/raw_reference.fa.gz.completed"),
#     params: url = lambda wildcards, output: genomes[wildcards.genome]["reference"]
#     wildcard_constraints:
#         genome = "((?!downloads).)*",
#     conda: "config/conda-download-env.yml"
#     threads: 1
#     shell: """
#         mkdir -p "$(dirname '{output.reference}')"
#         wget --progress=dot:mega -O '{output.reference}' '{params.url}'
#         """

# def reference_input_file_name(wildcards):
#     try:
#         genome_name = wildcards.genome
#         if corrected_genomes is not None and genome_name in corrected_genomes:
#             genome_name = corrected_genomes[genome_name]["source_genome"]
#         reference = genomes[genome_name]["reference"]

#         if reference.split('.')[-1] == "gz":
#             input_file_name = DATADIR + "downloads/" + genome_name + "/packed-reference/raw_reference.fa"
#         else:
#             input_file_name = DATADIR + "downloads/" + genome_name + "/reference/raw_reference.fa"
#         return input_file_name
#     except Exception:
#         traceback.print_exc()
#         sys.exit("Catched exception")

# localrules: link_reference
# rule link_reference:
#     input: file = reference_input_file_name,
#     output: file = GENOME_REFERENCE_FORMAT,
#     wildcard_constraints:
#         genome = "((?!downloads).)*",
#     threads: 1
#     shell: """
#         mkdir -p "$(dirname '{output.file}')"
#         ln -sr -T '{input.file}' '{output.file}'
#        """

#############################
###### Corrected Reads ######
#############################

# def correction_read_url_file_format(corrected_genome):
#     format = url_file_format(corrected_genomes[corrected_genome]["correction_short_reads"][0])

#     if format == "fasta":
#         return "fa"
#     elif corrected_genomes[corrected_genome].get("format", None) == "sra":
#         return "sra"
#     else:
#         return format

# localrules: download_correction_short_reads
# rule download_correction_short_reads:
#     output: file = DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.{format}",
#             completed = touch(DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.{format}.completed"),
#     params: url = lambda wildcards: corrected_genomes[wildcards.corrected_genome]["correction_short_reads"][int(wildcards.index)],
#             url_format = lambda wildcards: correction_read_url_file_format(wildcards.corrected_genome),
#             checksum = lambda wildcards: corrected_genomes[wildcards.corrected_genome]["correction_short_reads_checksum"] if "correction_short_reads_checksum" in corrected_genomes[wildcards.corrected_genome] else "",
#     wildcard_constraints:
#         format = "(fa|bam|sra)",
#         index = "\d+",
#         corrected_genome = "((?!corrected_reads).)*",
#     conda: "config/conda-download-env.yml"
#     threads: 1
#     shell: """
#         mkdir -p "$(dirname '{output.file}')"

#         if [ '{params.url_format}' != '{wildcards.format}' ]; then
#             echo "Error: url format '{params.url_format}' does not match format '{wildcards.format}' given by rule wildcard!"
#             exit 1
#         fi

#         if [ -z "{params.checksum}" ]; then
#             wget --progress=dot:mega -O '{output.file}' '{params.url}'
#         else
#             wget --progress=dot:mega -O '{output.file}' '{params.url}'
#             echo "Checksum given, but not supported yet"
#             exit 1
#         fi
#         """

# rule convert_correction_short_reads:
#     input: file = lambda wildcards: DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads." + correction_read_url_file_format(wildcards.corrected_genome),
#     output: file = DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.converted.fa",
#             completed = touch(DATADIR + "corrected_reads/{corrected_genome}/reads-{index}/reads.converted.fa.completed"),
#     params: file_format = lambda wildcards: correction_read_url_file_format(wildcards.corrected_genome),
#     wildcard_constraints:
#         index = "\d+",
#         genome = "((?!corrected_reads).)*",
#     conda: "config/conda-convert-reads-env.yml"
#     threads: 1
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
#     shell: """
#         if [ '{params.file_format}' == 'bam' ]; then
#             samtools fasta '{input.file}' > '{output.file}'
#         elif [ '{params.file_format}' == 'sra' ]; then
#             fastq-dump --stdout --fasta default '{input.file}' > '{output.file}'
#         else
#             ln -sr -T '{input.file}' '{output.file}'
#         fi
#         """

# rule combine_correction_short_reads:
#     input: files = lambda wildcards: expand(DATADIR + "corrected_reads/{{corrected_genome}}/reads-{index}/reads.converted.fa", index=range(len(corrected_genomes[wildcards.corrected_genome]["correction_short_reads"]))),
#     output: reads = CORRECTION_SHORT_READS_FORMAT,
#             completed = touch(CORRECTION_SHORT_READS_FORMAT + ".completed"),
#     params: input_list = lambda wildcards, input: "'" + "' '".join(input.files) + "'",
#     wildcard_constraints:
#         genome = "((?!corrected_reads).)*",
#     threads: 1
#     resources:
#         time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 360),
#         queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 360),
#     shell: "cat {params.input_list} > '{output.reads}'"

# rule ratatosk:
#     input:  correction_short_reads = CORRECTION_SHORT_READS_FORMAT,
#             long_reads = lambda wildcards: GENOME_READS_FORMAT.format(genome = corrected_genomes[wildcards.corrected_genome]["source_genome"]),
#             binary = RATATOSK_BINARY,
#     output: corrected_long_reads = DATADIR + "corrected_reads/{corrected_genome}/ratatosk/corrected_reads.fa",
#             completed = touch(DATADIR + "corrected_reads/{corrected_genome}/ratatosk/corrected_reads.fa.completed"),
#     wildcard_constraints:
#         corrected_genome = "((?!/ratatosk).)*",
#     threads: MAX_THREADS
#     resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 80000),
#                cpus = MAX_THREADS,
#                time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 2520),
#                queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 2520, 80000),
#                mail_type = "END",
#     shell: """
#         {input.binary} -v -c {threads} -s {input.correction_short_reads} -l {input.long_reads} -o {output.corrected_long_reads}
#         """

# localrules: select_read_corrector
# rule select_read_corrector:
#     input:  corrected_reads_from_read_corrector = DATADIR + "corrected_reads/{corrected_genome}/ratatosk/corrected_reads.fa",
#     output: corrected_reads = DATADIR + "corrected_reads/{corrected_genome}/corrected_reads.fa",
#     wildcard_constraints:
#         corrected_genome = "((?!/ratatosk).)*",
#     threads: 1
#     shell: "ln -sr '{input.corrected_reads_from_read_corrector}' '{output.corrected_reads}'"

#####################################
###### Homopolymer compression ######
#####################################

rule homopolymer_compress_reads:
    input:  reads = safe_format(GENOME_READS, homopolymer_compression = "none"),
            binary = HOMOPOLYMER_COMPRESS_RS_BINARY,
    output: reads = GENOME_READS,
    log:    HOCO_READS_LOG,
    params: threads = str(MAX_THREADS - 2),
    wildcard_constraints:
            homopolymer_compression = "yes",
    conda:  "config/conda-biopython-env.yml"
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 600),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 600, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 600, 1_000),
    # Use two threads, since the program uses two extra threads for IO.
    # More than two threads will likely not yield any speedup, as the application then becomes IO-bound.
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.binary}' --threads {params.threads} '{input.reads}' '{output.reads}' 2>&1 | tee '{log}'"

rule homopolymer_compress_reference:
    input:  reference = safe_format(GENOME_REFERENCE, homopolymer_compression = "none"),
            binary = HOMOPOLYMER_COMPRESS_RS_BINARY,
    output: reference = GENOME_REFERENCE,
    log:    HOCO_REFERENCE_LOG,
    params: threads = str(MAX_THREADS - 2),
    wildcard_constraints:
            homopolymer_compression = "yes",
    conda:  "config/conda-biopython-env.yml"
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 600),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 600, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 600, 1_000),
    # Use two threads, since the program uses two extra threads for IO.
    # More than two threads will likely not yield any speedup, as the application then becomes IO-bound.
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.binary}' --threads {params.threads} '{input.reference}' '{output.reference}' 2>&1 | tee '{log}'"

#############################
###### Read simulation ######
#############################

# Since fastk stores some outputs in the input dir, we need to link the input to the output dir.
# Otherwise parallel cluster executions with the same input dir will lock-conflict on Lustre.
localrules: link_fastk_reads
rule link_fastk_reads:
    input:  reads = GENOME_READS,
    output: reads = FASTK_INPUT_READS,
    shell: "ln -sr -T '{input.reads}' '{output.reads}'"

rule fastk_build_table_and_prof:
    input:  reads = FASTK_INPUT_READS,
            binary = FASTK_BINARY,
    output: table = FASTK_TABLE,
            profile = FASTK_PROFILE,
            hist = FASTK_HIST,
    params: tmp_dir = lambda wildcards, output: os.path.dirname(output.table)
    log:    os.path.join(FASTK_OUTPUT_DIR, "build_table.log")
    wildcard_constraints:
        read_source = "real",
        read_simulation_model_source = "none",
        uniquify_ids = "no",
        fastk_k = "(([4-9][0-9])|(1[0-1][0-9])|(12[0-8]))", # >= 40, <= 128
    threads: MAX_THREADS,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 120, 1_000),
    shell:  """
        '{input.binary}' -v -T{threads} -P'{params.tmp_dir}' -t1 -p -k{wildcards.fastk_k} '{input.reads}' 2>&1 | tee '{log}'
        """

rule histex_evaluation_table:
    input:  hist = FASTK_HIST,
            binary = FASTK_HISTEX_BINARY,
    log:    histogram = FASTK_HISTEX_EVALUATION,
    wildcard_constraints:
        read_source = "real",
        read_simulation_model_source = "none",
        uniquify_ids = "no",
        fastk_k = "(([4-9][0-9])|(1[0-1][0-9])|(12[0-8]))", # >= 40, <= 128
    threads: 1,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 10),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 10, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 10, 1_000),
    shell:  """
        '{input.binary}' -h10000 '{input.hist}' 2>&1 | tee '{log.histogram}'
        '{input.binary}' -k -h10000 '{input.hist}' 2>&1 | tee -a '{log.histogram}'
        """


rule symmex_table:
    input:  table = FASTK_TABLE,
            binary = FASTK_SYMMEX_BINARY,
    output: table = FASTK_SYMMETRIC_TABLE,
    params: tmp_dir = lambda wildcards, output: os.path.dirname(output.table)
    log:    os.path.join(FASTK_OUTPUT_DIR, "symmex_table.log")
    wildcard_constraints:
        read_source = "real",
        read_simulation_model_source = "none",
        uniquify_ids = "no",
        fastk_k = "(([4-9][0-9])|(1[0-1][0-9])|(12[0-8]))", # >= 40, <= 128
    threads: MAX_THREADS,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 1_000),
    shell:  """
        '{input.binary}' -v -T{threads} -P'{params.tmp_dir}' '{input.table}' '{output.table}' 2>&1 | tee '{log}'
        """

# Since himodel stores some outputs in the input dir, we need to link the input to the output dir.
# Otherwise parallel cluster executions with the same input dir will lock-conflict on Lustre.
localrules: link_himodel_symmetric_table
rule link_himodel_symmetric_table:
    input:  table = FASTK_SYMMETRIC_TABLE,
    output: table = HIMODEL_INPUT_SYMMETRIC_TABLE,
    params: input_filename = lambda wildcards, input: os.path.basename(input.table),
            input_dirname = lambda wildcards, input: os.path.dirname(input.table),
            output_filename = lambda wildcards, output: os.path.basename(output.table),
            output_dirname = lambda wildcards, output: os.path.dirname(output.table),
    shell:  """
        ln -sr -T '{input.table}' '{output.table}'
        
        for INPUT_FILE_NAME in $(ls '{params.input_dirname}/.{params.input_filename}.'* | xargs -n 1 basename); do
            INPUT_FILE='{params.input_dirname}'/"${{INPUT_FILE_NAME}}"
            OUTPUT_FILE_NAME=${{INPUT_FILE_NAME/{params.input_filename}/{params.output_filename}}}
            OUTPUT_FILE='{params.output_dirname}'/"${{OUTPUT_FILE_NAME}}"
            ln -sr -T "${{INPUT_FILE}}" "${{OUTPUT_FILE}}"
        done
        """

# Since himodel stores some outputs in the input dir, we need to link the input to the output dir.
# Otherwise parallel cluster executions with the same input dir will lock-conflict on Lustre.
localrules: link_himodel_profile
rule link_himodel_profile:
    input:  profile = FASTK_PROFILE,
    output: profile = HIMODEL_INPUT_PROFILE,
    params: input_filename = lambda wildcards, input: os.path.basename(input.profile),
            input_filename_pidx = lambda wildcards, input: os.path.basename(input.profile)[:-4] + "pidx",
            input_dirname = lambda wildcards, input: os.path.dirname(input.profile),
            output_filename = lambda wildcards, output: os.path.basename(output.profile),
            output_filename_pidx = lambda wildcards, output: os.path.basename(output.profile)[:-4] + "pidx",
            output_dirname = lambda wildcards, output: os.path.dirname(output.profile),
    shell:  """
        ln -sr -T '{input.profile}' '{output.profile}'
        
        for INPUT_FILE_NAME in $(ls '{params.input_dirname}/.{params.input_filename}.'* | xargs -n 1 basename); do
            INPUT_FILE='{params.input_dirname}'/"${{INPUT_FILE_NAME}}"
            OUTPUT_FILE_NAME=${{INPUT_FILE_NAME/{params.input_filename}/{params.output_filename}}}
            OUTPUT_FILE='{params.output_dirname}'/"${{OUTPUT_FILE_NAME}}"
            ln -sr -T "${{INPUT_FILE}}" "${{OUTPUT_FILE}}"
        done
        for INPUT_FILE_NAME in $(ls '{params.input_dirname}/.{params.input_filename_pidx}.'* | xargs -n 1 basename); do
            INPUT_FILE='{params.input_dirname}'/"${{INPUT_FILE_NAME}}"
            OUTPUT_FILE_NAME=${{INPUT_FILE_NAME/{params.input_filename_pidx}/{params.output_filename_pidx}}}
            OUTPUT_FILE='{params.output_dirname}'/"${{OUTPUT_FILE_NAME}}"
            ln -sr -T "${{INPUT_FILE}}" "${{OUTPUT_FILE}}"
        done
        """

def hisim_model_input_prefix(wildcards, input):
    try:
        result = input.symmetric_table[:-5]
        if result != input.profile[:-5]:
            raise Exception(f"mismatching symmetric table and profile files (need to match except for the extension): '{input.symmetric_table}' and '{input.profile}'")
        return os.path.basename(result)
    except:
        traceback.print_exc()
        sys.exit("Catched exception")

rule hisim_model:
    input:  symmetric_table = HIMODEL_INPUT_SYMMETRIC_TABLE,
            profile = HIMODEL_INPUT_PROFILE,
            binary = HISIM_MODEL_BINARY,
    output: model = HIMODEL_MODEL,
    log:    os.path.join(HIMODEL_OUTPUT_DIR, "himodel.log")
    params: working_directory = HIMODEL_OUTPUT_DIR,
            input_binary = lambda wildcards, input: os.path.abspath(input.binary),
            input_prefix = hisim_model_input_prefix,
            log = "himodel.log",
    wildcard_constraints:
        homopolymer_compression = "none",
        read_source = "real",
        read_simulation_model_source = "none",
        read_downsampling_factor = "none",
        uniquify_ids = "no",
        fastk_k = "(([4-9][0-9])|(1[0-1][0-9])|(12[0-8]))", # >= 40, <= 128
        himodel_kmer_threshold = "[1-9][0-9]*",
        himodel_min_valid = "[1-9][0-9]*",
        himodel_max_valid = "[1-9][0-9]*",
    threads: MAX_THREADS,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = MAX_THREADS,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 1_000),
    shell:  """
        cd '{params.working_directory}'
        '{params.input_binary}' -v -T{threads} -g{wildcards.himodel_min_valid}:{wildcards.himodel_max_valid} -e{wildcards.himodel_kmer_threshold} '{params.input_prefix}' 2>&1 | tee '{params.log}'
        """

def hisim_sim_genome_reference(wildcards):
    try:
        read_simulation_model_source = wildcards.read_simulation_model_source
        if read_simulation_model_source == "none":
            read_simulation_model_source = wildcards.genome

        genome = genomes[read_simulation_model_source]
        genome_arguments = genome.setdefault("genome_arguments", {})
        
        filter_nw = genome_arguments.setdefault("filter_nw", "no")
        retain_cm = genome_arguments.setdefault("retain_cm", "no")
        filter_plasmids = genome_arguments.setdefault("filter_plasmids", "no")
        return safe_format(GENOME_REFERENCE, genome = read_simulation_model_source, filter_nw = filter_nw, retain_cm = retain_cm, filter_plasmids = filter_plasmids)
    except:
        traceback.print_exc()
        sys.exit("Catched exception")

def hisim_sim_model(wildcards):
    try:
        read_source = wildcards.read_source
        read_simulation_model_source = wildcards.read_simulation_model_source
        if read_simulation_model_source == "none":
            read_simulation_model_source = wildcards.genome

        genome = genomes[read_simulation_model_source]
        genome_arguments = genome.setdefault("genome_arguments", {})
        assert "fastk_k" in genome_arguments, f"genome_arguments of '{genome}' misses 'fastk_k'"
        assert "himodel_kmer_threshold" in genome_arguments, f"genome_arguments of '{genome}' misses 'himodel_kmer_threshold'"
        assert "himodel_min_valid" in genome_arguments, f"genome_arguments of '{genome}' misses 'himodel_min_valid'"
        assert "himodel_max_valid_low" in genome_arguments, f"genome_arguments of '{genome}' misses 'himodel_max_valid_low'"
        assert "himodel_max_valid_medium" in genome_arguments, f"genome_arguments of '{genome}' misses 'himodel_max_valid_medium'"
        assert "himodel_max_valid_high" in genome_arguments, f"genome_arguments of '{genome}' misses 'himodel_max_valid_high'"

        fastk_k = genome_arguments["fastk_k"]
        himodel_kmer_threshold = genome_arguments["himodel_kmer_threshold"]
        himodel_min_valid = genome_arguments["himodel_min_valid"]

        if read_source.endswith("_low"):
            himodel_max_valid = genome_arguments["himodel_max_valid_low"]
        elif read_source.endswith("_medium"):
            himodel_max_valid = genome_arguments["himodel_max_valid_medium"]
        elif read_source.endswith("_high"):
            himodel_max_valid = genome_arguments["himodel_max_valid_high"]
        else:
            raise Exception(f"read_source does not end with an identifier for the upper bound to use '{read_source}'")

        return safe_format(HIMODEL_MODEL, genome = read_simulation_model_source, read_source = "real", read_simulation_model_source = "none", read_downsampling_factor = "none", uniquify_ids = "no", fastk_k = fastk_k, himodel_kmer_threshold = himodel_kmer_threshold, himodel_min_valid = himodel_min_valid, himodel_max_valid = himodel_max_valid)
    except:
        traceback.print_exc()
        sys.exit("Catched exception")

def hisim_sim_params(wildcards):
    try:
        clean_read_source = wildcards.read_source
        if clean_read_source.endswith("_low"):
            clean_read_source = clean_read_source[:-4]
        elif clean_read_source.endswith("_medium"):
            clean_read_source = clean_read_source[:-7]
        elif clean_read_source.endswith("_high"):
            clean_read_source = clean_read_source[:-5]

        if clean_read_source == "hisim_test":
            return "-c40.0"
        elif clean_read_source == "hisim_human":
            return "-c25.0"
        elif clean_read_source == "hisim_human_haploid":
            return "-c25.0"
        elif clean_read_source == "hisim_0_7":
            return ""
        elif clean_read_source == "hisim_haploid":
            return ""
        elif clean_read_source == "hisim_test_x0":
            return "-x0"
        elif clean_read_source == "hisim_test_x1000":
            return "-x1000"
        elif clean_read_source == "hisim_test_x2000":
            return "-x2000"
        elif clean_read_source == "hisim_test_x3000":
            return "-x3000"
        elif clean_read_source == "hisim_test_x4000":
            return "-x4000"
        else:
            raise Exception(f"Unknown read source: {wildcards.read_source} (clean: {clean_read_source})")
    except:
        traceback.print_exc()
        sys.exit("Catched exception")

def hisim_ploidy_tree(wildcards):
    try:
        clean_read_source = wildcards.read_source
        if clean_read_source.endswith("_low"):
            clean_read_source = clean_read_source[:-4]
        elif clean_read_source.endswith("_medium"):
            clean_read_source = clean_read_source[:-7]
        elif clean_read_source.endswith("_high"):
            clean_read_source = clean_read_source[:-5]

        if clean_read_source == "hisim_test":
            return "0.0,0.2"
        elif clean_read_source == "hisim_human":
            return "0.0,0.66"
        elif clean_read_source == "hisim_human_haploid":
            return "0.0"
        elif clean_read_source == "hisim_0_7":
            return "0.0,0.7"
        elif clean_read_source == "hisim_haploid":
            return "0.0"
        elif clean_read_source.startswith("hisim_test_x"):
            return "0.0"
        else:
            raise Exception(f"Unknown read source: {wildcards.read_source} (clean: {clean_read_source})")
    except:
        traceback.print_exc()
        sys.exit("Catched exception")

rule hisim_sim:
    input:  reference = hisim_sim_genome_reference,
            model = hisim_sim_model,
            binary = HISIM_SIM_BINARY,
    output: reads = GENOME_READS,
            haplotype = safe_format(HISIM_HAPLOTYPE, haplotype_index = "1"),
    params: sim_params = hisim_sim_params,
            output_prefix = lambda wildcards, output: os.path.join(os.path.dirname(output.reads), os.path.basename(output.reads).replace(".fa", "")),
            ploidy_tree = hisim_ploidy_tree,
    log:    os.path.join(GENOME_READS_DIR, "hisim.log"),
    wildcard_constraints:
        read_source = "hisim_.*",
        read_downsampling_factor = "none",
        uniquify_ids = "no",
    threads: 1,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 1_000),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 1_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 1_000),
    shell:  """
        '{input.binary}' -v '{input.reference}' '{input.model}' -o'{params.output_prefix}' {params.sim_params} -p{params.ploidy_tree} -fh -r3541529 -U 2>&1 | tee '{log}'
        ln -sr -T '{params.output_prefix}.fasta' '{output.reads}'
        """

###############################
###### Read downsampling ######
###############################

rule downsample_reads:
    input:  reads = safe_format(GENOME_READS, read_downsampling_factor = "none"),
            script = DOWNSAMPLE_FASTA_READS_SCRIPT,
    output: reads = GENOME_READS,
    wildcard_constraints:
            read_downsampling_factor = "0.[0-9]+",
            homopolymer_compression = "none",
    conda:  "config/conda-biopython-env.yml"
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.script}' '{input.reads}' '{output.reads}' {wildcards.read_downsampling_factor}"

###########################################
###### Convert to single-lined fasta ######
###########################################

rule convert_reads_to_single_lined_fasta:
    input:  reads = GENOME_READS,
    output: reads = GENOME_SINGLE_LINE_READS,
    conda:  "config/conda-seqtk-env.yml"
    shell:  "${{CONDA_PREFIX}}/bin/time -v seqtk seq -AU '{input.reads}' > '{output.reads}'"

############################################
###### Filter NW entries in reference ######
############################################

rule filter_nw:
    input:  reference = safe_format(GENOME_REFERENCE, filter_nw = "no", filter_plasmids = "no"),
            script = FILTER_NW_FROM_REFERENCE_SCRIPT,
    output: reference = GENOME_REFERENCE,
    wildcard_constraints:
            homopolymer_compression = "none",
            filter_nw = "yes",
            retain_cm = "no",
    conda:  "config/conda-biopython-env.yml"
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.script}' '{input.reference}' '{output.reference}'"

#################################################
###### Retain only CM entries in reference ######
#################################################

rule retain_cm:
    input:  reference = safe_format(GENOME_REFERENCE, retain_cm = "no", filter_plasmids = "no"),
            script = RETAIN_CM_IN_REFERENCE_SCRIPT,
    output: reference = GENOME_REFERENCE,
    wildcard_constraints:
            homopolymer_compression = "none",
            retain_cm = "yes",
            filter_nw = "no",
    conda:  "config/conda-biopython-env.yml"
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.script}' '{input.reference}' '{output.reference}'"

#################################################
###### Filter plasmid entries in reference ######
#################################################

rule filter_plasmids:
    input:  reference = safe_format(GENOME_REFERENCE, retain_cm = "no", filter_nw = "no"),
            script = FILTER_PLASMIDS_IN_REFERENCE_SCRIPT,
    output: reference = GENOME_REFERENCE,
    wildcard_constraints:
            homopolymer_compression = "none",
            retain_cm = "no",
            filter_nw = "no",
            filter_plasmids = "yes",
    conda:  "config/conda-biopython-env.yml"
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.script}' '{input.reference}' '{output.reference}'"

######################################
###### Compute reference length ######
######################################

rule compute_reference_length:
    input:  reference = GENOME_REFERENCE,
            script = COMPUTE_GENOME_REFERENCE_LENGTH_SCRIPT,
    output: reference_length = GENOME_REFERENCE_LENGTH,
    conda:  "config/conda-biopython-env.yml"
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.script}' '{input.reference}' '{output.reference_length}'"

######################################
###### Input Genome Preparation ######
######################################

rule separate_linear_and_circular:
    input: filtered = DATADIR + "{genome}/filtered.fna", verified = DATADIR + "{genome}/is_genome_verified.log", binary = RUST_BINARY,
    output: circular = DATADIR + "{genome}/circular.fna", linear = DATADIR + "{genome}/linear.fna", log = DATADIR + "{genome}/separate_linear_and_circular.log"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "cp '{input.filtered}' '{output.linear}'; data/target/release/cli circularise-genome --input '{input.filtered}' 2>&1 --output '{output.circular}' | tee '{output.log}'"

rule verify_genome:
    input: file = DATADIR + "{dir}/filtered.fna", binary = RUST_BINARY,
    output: log = DATADIR + "{dir}/is_genome_verified.log"
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: DATADIR + "target/release/cli verify-genome --input '{input.file}' 2>&1 | tee '{output.log}'"

rule filter_genome:
    input: file = DATADIR + "{dir}/raw.fna", binary = RUST_BINARY,
    output: file = DATADIR + "{dir}/filtered.fna", genome_name = DATADIR + "{dir}/name.txt", log = DATADIR + "{dir}/filtered.log"
    params: retain = lambda wildcards: "--retain '" + experiments_bcalm2[wildcards.dir]["filter_retain"] + "'" if "filter_retain" in experiments_bcalm2[wildcards.dir] else ""
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: DATADIR + "target/release/cli filter --input '{input.file}' --output '{output.file}' --extract-name '{output.genome_name}' {params.retain} 2>&1 | tee '{output.log}'"

###################
###### QUAST ######
###################

def get_quast_extra_arguments_from_wildcards(wildcards):
    try:
        if wildcards.quast_mode == "hicanu":
            return ""
        elif wildcards.quast_mode == "hicanu_alignments":
            return "--skip-unaligned-mis-contigs --min-alignment 10000 --min-identity 98.0 --extensive-mis-size 5000 --min-contig 50000"
        elif wildcards.quast_mode == "hicanu_misassemblies":
            return "--min-alignment 20000 --extensive-mis-size 500000 --min-identity 90"
        elif wildcards.quast_mode == "hicanu_misassemblies_strict":
            return "--min-alignment 20000 --extensive-mis-size 100000 --min-identity 90"
        elif wildcards.quast_mode == "normal":
            return "--fragmented"
        elif wildcards.quast_mode == "hicanu_hoco":
            return "--minimap-hoco"
        elif wildcards.quast_mode == "hicanu_alignments_hoco":
            return "--skip-unaligned-mis-contigs --min-alignment 10000 --min-identity 98.0 --extensive-mis-size 5000 --min-contig 50000 --minimap-hoco"
        elif wildcards.quast_mode == "hicanu_misassemblies_hoco":
            return "--min-alignment 20000 --extensive-mis-size 500000 --min-identity 90 --minimap-hoco"
        elif wildcards.quast_mode == "hicanu_misassemblies_strict_hoco":
            return "--min-alignment 20000 --extensive-mis-size 100000 --min-identity 90 --minimap-hoco"
        elif wildcards.quast_mode == "normal_hoco":
            return "--fragmented --minimap-hoco"
        elif wildcards.quast_mode == "hicanu_wrapped_hoco":
            return "--minimap-hoco-wrapped"
        elif wildcards.quast_mode == "hicanu_alignments_wrapped_hoco":
            return "--skip-unaligned-mis-contigs --min-alignment 10000 --min-identity 98.0 --extensive-mis-size 5000 --min-contig 50000 --minimap-hoco-wrapped"
        elif wildcards.quast_mode == "hicanu_misassemblies_wrapped_hoco":
            return "--min-alignment 20000 --extensive-mis-size 500000 --min-identity 90 --minimap-hoco-wrapped"
        elif wildcards.quast_mode == "hicanu_misassemblies_strict_wrapped_hoco":
            return "--min-alignment 20000 --extensive-mis-size 100000 --min-identity 90 --minimap-hoco-wrapped"
        elif wildcards.quast_mode == "normal_wrapped_hoco":
            return "--fragmented --minimap-hoco-wrapped"
        else:
            raise Exception(f"Unknown quast_mode: {wildcards.quast_mode}")
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_quast_references_from_wildcards(wildcards):
    try:
        if wildcards.read_source == "real":
            return "-r '" + GENOME_REFERENCE.format(**wildcards) + "'"
        else:
            ploidy_tree = hisim_ploidy_tree(wildcards)
            ploidy_count = len(ploidy_tree.split(","))

            files = []
            for haplotype_index in range(1, ploidy_count + 1):
                files.append(safe_format(HISIM_HAPLOTYPE, haplotype_index = haplotype_index, uniquify_ids = "no").format(**wildcards))

            return "-r '" + "' -r '".join(files) + "'"
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule run_quast:
    input:  contigs = ASSEMBLED_CONTIGS,
            reference = GENOME_REFERENCE,
            hisim_haplotype = lambda wildcards: [] if wildcards.read_source == "real" else safe_format(HISIM_HAPLOTYPE, haplotype_index = "1", uniquify_ids = "no"),
            script = QUAST_BINARY,
            homopolymer_compress_rs_binary = HOMOPOLYMER_COMPRESS_RS_BINARY,
            minimap2_homopolymer_decompression_binary = MINIMAP2_HOMOPOLYMER_DECOMPRESSION_BINARY,
    output: directory = directory(QUAST_OUTPUT_DIR),
            eaxmax_csv = os.path.join(QUAST_OUTPUT_DIR, "aligned_stats/EAxmax_plot.csv"),
    log:    minimap = os.path.join(QUAST_OUTPUT_DIR, "contigs_reports", "contigs_report_contigs.stderr")
    params: extra_arguments = get_quast_extra_arguments_from_wildcards,
            references = get_quast_references_from_wildcards,
    conda: "config/conda-quast-env.yml"
    threads: 14,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 50_000),
               cpus = 14,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 120),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 120, 50_000),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 120, 50_000),
    shell: "${{CONDA_PREFIX}}/bin/time -v {input.script} {params.extra_arguments} -t {threads} --hoco-binary '{input.homopolymer_compress_rs_binary}' --hodeco-binary '{input.minimap2_homopolymer_decompression_binary}' --no-html --large -o '{output.directory}' {params.references} '{input.contigs}'"

##################################
###### Resources evaluation ######
##################################

def get_evaluate_resources_inputs(wildcards):
    try:
        result = {}
        assembler_arguments = wildcards.assembler_arguments
        assembler_arguments = parse.parse(ASSEMBLER_ARGUMENT_STRINGS[wildcards.assembler], wildcards.assembler_arguments).named

        if wildcards.assembler == "wtdbg2":
            # hoco resources
            if assembler_arguments["hodeco_consensus"] != "none":
                assembly_hoco = "yes"
                result["hoco"] = safe_format(HOCO_READS_LOG, homopolymer_compression = "yes").format(**wildcards)
                result["hodeco"] = safe_format(os.path.join(WTDBG2_OUTPUT_DIR, "transform_ctg_lay_hodeco_simple.log"), **assembler_arguments).format(**wildcards)
            else:
                assembly_hoco = wildcards.homopolymer_compression

            # injections
            has_injections = False
            if assembler_arguments["tig_injection"] != "none":
                result["trivial_omnitigs"] = safe_format(safe_format(os.path.join(WTDBG2_INJECTABLE_CONTIG_DIR, "compute_injectable_contigs.log"), homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)
                has_injections = True
            elif assembler_arguments["frg_injection"] != "none":
                result["trivial_omnitigs"] = safe_format(safe_format(os.path.join(WTDBG2_INJECTABLE_FRAGMENT_CONTIG_DIR, "compute_injectable_contigs.log"), homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)
                has_injections = True

            if has_injections:
                # two-step assembly
                result["assembly"] = safe_format(safe_format(WTDBG2_LOG, tig_injection = "none", frg_injection = "none", frg_injection_stage = "none", homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)
                result["wtdbg2_extract"] = safe_format(safe_format(WTDBG2_EXTRACT_LOG, subfile = "ctg.lay", tig_injection = "none", frg_injection = "none", frg_injection_stage = "none", homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)
                result["contig_assembly"] = safe_format(safe_format(WTDBG2_LOG, homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)
            else:
                # one-step assembly
                result["assembly"] = safe_format(safe_format(WTDBG2_LOG, homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)
                # extract
                result["wtdbg2_extract"] = safe_format(safe_format(WTDBG2_EXTRACT_LOG, subfile = "ctg.lay", homopolymer_compression = assembly_hoco, hodeco_consensus = "none"), **assembler_arguments).format(**wildcards)


            # consensus
            result["wtdbg2_consensus"] = safe_format(WTDBG2_CONSENSUS_LOG, **assembler_arguments).format(**wildcards)
        elif wildcards.assembler == "hifiasm":
            # assembly
            result["assembly"] = safe_format(safe_format(HIFIASM_LOG, contig_algorithm = "builtin")).format(**wildcards)

            # injections
            if assembler_arguments["contig_algorithm"].startswith("trivial_omnitigs_"):
                result["trivial_omnitigs"] = safe_format(os.path.join(HIFIASM_OUTPUT_DIR, "compute_injectable_contigs.log"), **assembler_arguments).format(**wildcards)
        else:
            # assembly
            result["assembly"] = ASSEMBLY_LOG.format(**wildcards)

        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def decode_time(string):
    try:
        if string.count(':') == 2:
            hours, minutes, seconds = map(float, string.split(':'))
            time = datetime.timedelta(hours = hours, minutes = minutes, seconds = seconds)
        elif string.count(':') == 1:
            minutes, seconds = map(float, string.split(':'))
            time = datetime.timedelta(minutes = minutes, seconds = seconds)
        else:
            raise Exception(f"unknown time string '{string}'")
        return time.total_seconds()
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

rule evaluate_resources:
    input:  files = lambda wildcards: get_evaluate_resources_inputs(wildcards).values(),
    output: file = RESOURCES_EVALUATION,
    params: file_map = get_evaluate_resources_inputs,
    threads: 1,
    run:
        result = {}
        for key, input_file_name in params.file_map.items():
            with open(input_file_name, 'r') as input_file:
                values = {"time": 0, "mem": 0}
                for line in input_file:
                    if "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                        line = line.replace("Elapsed (wall clock) time (h:mm:ss or m:ss):", "").strip()
                        values["time"] = decode_time(line) + values["time"]
                    elif "Maximum resident set size" in line:
                        values["mem"] = max(int(line.split(':')[1].strip()), values["mem"])

                assert "time" in values, f"No time found in {input_file_name}"
                assert "mem" in values, f"No mem found in {input_file_name}"
                result[key] = values

        sum_time = sum([values["time"] for values in result.values()])
        max_mem = max([values["mem"] for values in result.values()])
        result["total"] = {
            "time": sum_time,
            "mem": max_mem,
        }

        with open(output.file, 'w') as output_file:
            json.dump(result, output_file)


#############################
###### ContigValidator ######
#############################

rule run_contig_validator:
    input:  contig_validator_dir = CONTIG_VALIDATOR_DIR,
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
        cd '{input.contig_validator_dir}'
        # The abundance-min here has nothing to do with the abundance_min from bcalm2
        bash run.sh -suffixsave 0 -abundance-min 1 -kmer-size {wildcards.k} -r '../../{input.reference}' -a '../../{output.result}' -i '../../{input.reads}'
        """

#####################
###### Bandage ######
#####################

rule convert_bcalm2_output_to_gfa:
    input:  fa = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.fa",
            converter = CONVERT_TO_GFA_BINARY,
    output: gfa = DATADIR + "{dir}/{file}.k{k}-a{abundance_min}.bcalm2.gfa",
    threads: 1
    shell:  "${{CONDA_PREFIX}}/bin/time -v '{input.converter}' {input.fa} {output.gfa} {wildcards.k}"

rule bandage:
    input: "{file}.gfa"
    output: "{file}.bandage.png"
    conda: "config/conda-bandage-env.yml"
    threads: 1
    shell: "Bandage image {input} {output} --width 1000 --height 1000"

#================================================
#=== Downloads ==================================
#================================================

def escape_dirname(raw):
    try:
        assert type(raw) is str, "type of raw dirname is not str"
        assert "¤" not in raw, "raw dirname contains ¤"
        assert raw != "None", "raw dirname is textual None"

        escaped = raw.replace("/", "¤a").replace("?", "¤b").replace("=", "¤c").replace("&", "¤d")

        result = ""
        while len(escaped) > 200:
            result += escaped[:200]
            escaped = escaped[200:]
            if len(escaped) > 0:
                result += "/"
        result += escaped

        assert result is not None, "escaped dirname is None"
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def unescape_dirname(escaped):
    try:
        assert type(escaped) is str
        assert escaped != "None", "escaped is textual None"

        raw = escaped.replace("/", "").replace("¤a", "/").replace("¤b", "?").replace("¤c", "=").replace("¤d", "&")
        assert raw is not None
        return raw
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def checksum_url(url):
    try:
        assert type(url) is str
        assert url != "None", "url is textual None"

        if "hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/" in url:
            return "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt"

        url = urlparse(url)
        dirname = os.path.dirname(url.path)
        url = url._replace(path = os.path.join(dirname, "md5checksums.txt"))
        url = url.geturl()
        assert url is not None
        return url
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

# download files

localrules: download_fa_file
rule download_fa_file:
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fa"),
            checksum_file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "md5checksums.txt"),
    params: url = lambda wildcards: unescape_dirname(wildcards.url),
            checksum_url = lambda wildcards: checksum_url(unescape_dirname(wildcards.url)),
    wildcard_constraints:
            url = "http.*\./?(f/?a|f/?n/?a|f/?a/?s/?t/?a)"
    shell:  """
        wget --progress=dot:mega -O '{output.file}' '{params.url}'
        wget --progress=dot:mega -O '{output.checksum_file}' '{params.checksum_url}'

        CHECKSUM=$(md5sum '{output.file}' | cut -f1 -d' ' | sed 's/[\]//g')
        cat '{output.checksum_file}' | grep "$CHECKSUM"
    """

#use rule download_fa_file as download_fa_gz_file with:
#    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fa.gz"),
#    wildcard_constraints:
#            url = ".*\./?(f/?a|f/?n/?a|f/?a/?s/?t/?a)/?\./?g/?z"

localrules: download_fa_gz_file
rule download_fa_gz_file:
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fa.gz"),
            checksum_file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "md5checksums.txt"),
    params: url = lambda wildcards: unescape_dirname(wildcards.url),
            checksum_url = lambda wildcards: checksum_url(unescape_dirname(wildcards.url)),
    wildcard_constraints:
            url = "http.*\./?((f/?a)|(f/?n/?a)|(f/?a/?s/?t/?a))/?\./?g/?z"
    shell:  """
        wget --progress=dot:mega -O '{output.file}' '{params.url}'
        wget --progress=dot:mega -O '{output.checksum_file}' '{params.checksum_url}'

        CHECKSUM=$(md5sum '{output.file}' | cut -f1 -d' ' | sed 's/[\]//g')
        echo $CHECKSUM
        cat '{output.checksum_file}' | grep "$CHECKSUM"
    """

localrules: download_fastq_gz_file
rule download_fastq_gz_file:
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fastq.gz"),
            #checksum_file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "md5checksums.txt"),
    params: url = lambda wildcards: unescape_dirname(wildcards.url),
            #checksum_url = lambda wildcards: checksum_url(unescape_dirname(wildcards.url)),
    wildcard_constraints:
            url = "http.*\./?((f/?q)|(f/?n/?q)|(f/?a/?s/?t/?q))/?\./?g/?z"
    shell:  """
        wget --progress=dot:mega -O '{output.file}' '{params.url}'
    """

localrules: download_sra_file
rule download_sra_file:
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.sra"),
    params: url = lambda wildcards: unescape_dirname(wildcards.url),
    wildcard_constraints:
            url = "http.*(((S|D)/?R/?R)|((S|D)/?R/?A))[0-9/\.]+"
    shell:  """
        wget --progress=dot:mega -O '{output.file}' '{params.url}'
    """

# convert files

rule convert_fastq_download:
    input:  file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fastq"),
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fa"),
    conda:  "config/conda-convert-reads-env.yml"
    shell:  """
        bioawk -c fastx '{{ print ">" $name "\\n" $seq }}' '{input.file}' > '{output.file}'
    """

rule convert_sra_download:
    input:  file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.sra"),
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "file.fa"),
    conda:  "config/conda-convert-reads-env.yml"
    shell:  "fastq-dump --stdout --fasta default '{input.file}' > '{output.file}'"

rule extract_download:
    input:  file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "{file}.gz"),
    output: file = os.path.join(DOWNLOAD_ROOTDIR, "file", "{url}", "{file}"),
    params: working_directory = lambda wildcards, input: os.path.dirname(input.file),
    wildcard_constraints:
            file = "[^/]*(?<!\.gz)",
    conda:  "config/conda-extract-env.yml"
    shell:  "cd '{params.working_directory}'; gunzip -k {wildcards.file}.gz"

def get_genome_url(wildcards):
    try:
        genome = str(wildcards.genome)
        assert genome in genomes, f"Genome not found: {genome}"

        result = genomes[genome]["reference"]

        if result is None:
            raise Exception("Did not find genome url")
        else:
            assert result != "None", "result is textual None."
            return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

def get_genome_reads_urls(wildcards):
    try:
        genome = str(wildcards.genome)
        assert genome in genomes, f"Genome not found: {genome}"

        result = genomes[genome]["reads"]
        
        assert result != "None", "result is textual None."
        if result is None:
            raise Exception("Did not find genome url")
        
        if type(result) is str:
            result = [result]

        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")

# general linking rules

localrules: download_reference_genome
rule download_reference_genome:
    input:  file = lambda wildcards: os.path.join(DOWNLOAD_ROOTDIR, "file", escape_dirname(get_genome_url(wildcards)), "file.fa"),
    output: file = GENOME_REFERENCE,
    wildcard_constraints:
            homopolymer_compression = "none",
            filter_nw = "no",
            retain_cm = "no",
    shell: "ln -sr -T '{input.file}' '{output.file}'"

rule download_genome_reads:
    input:  files = lambda wildcards: [os.path.join(DOWNLOAD_ROOTDIR, "file", escape_dirname(url), "file.fa") for url in get_genome_reads_urls(wildcards)],
    output: file = GENOME_READS,
    wildcard_constraints:
            read_source = "real",
            read_simulation_model_source = "none",
            read_downsampling_factor = "none",
            homopolymer_compression = "none",
            uniquify_ids = "no",
    params: input_files = lambda wildcards, input: "'" + "' '".join(input.files) + "'"
    threads: 1,
    resources: mem_mb = lambda wildcards: compute_genome_mem_mb_from_wildcards(wildcards, 100),
               cpus = 1,
               time_min = lambda wildcards: compute_genome_time_min_from_wildcards(wildcards, 60),
               queue = lambda wildcards: compute_genome_queue_from_wildcards(wildcards, 60, 100),
               cluster = lambda wildcards: compute_genome_cluster_from_wildcards(wildcards, 60, 100),
    shell:  "cat {params.input_files} > '{output.file}'"

rule uniquify_ids:
    input:  reads = safe_format(GENOME_READS, uniquify_ids = "no"),
            script = UNIQUIFY_IDS_SCRIPT,
    output: reads = GENOME_READS,
    log:    log = UNIQUIFY_IDS_LOG,
    wildcard_constraints:
            read_downsampling_factor = "none",
            homopolymer_compression = "none",
            uniquify_ids = "yes",
    conda: "config/conda-uniquify-env.yml"
    threads: 1
    shell: "python3 '{input.script}' '{input.reads}' '{output.reads}' 2>&1 | tee '{log.log}'"

###########################
###### Installations ######
###########################

localrules: fetch_rust
rule fetch_rust:
    input:  sources = RUST_SOURCES,
    output: is_rust_fetched_marker = touch(IS_RUST_FETCHED_MARKER),
    log:    log = os.path.join(RUST_DIR, "fetch.log"),
    conda: "config/conda-rust-env.yml"
    threads: 1
    shell: "cargo fetch --manifest-path 'implementation/Cargo.toml' 2>&1 | tee '{log.log}'"

rule test_rust:
    input:  is_rust_fetched_marker = IS_RUST_FETCHED_MARKER,
    output: is_rust_tested_marker = touch(IS_RUST_TESTED_MARKER),
    log:    log = os.path.join(RUST_DIR, "test.log"),
    params: rust_dir = RUST_DIR,
    conda: "config/conda-rust-env.yml"
    threads: BUILD_THREADS
    resources: mem_mb = 4000,
               cpus = BUILD_THREADS,
               time_min = 30,
    shell: "cargo test -j {threads} --target-dir '{params.rust_dir}' --manifest-path 'implementation/Cargo.toml' --offline 2>&1 | tee '{log.log}'"

rule build_rust_release:
    input:  is_rust_tested_marker = IS_RUST_TESTED_MARKER,
    output: binary = RUST_BINARY,
    log:    log = os.path.join(RUST_DIR, "build_release.log"),
    params: rust_dir = RUST_DIR,
    conda: "config/conda-rust-env.yml"
    threads: BUILD_THREADS
    resources: mem_mb = 4000,
               cpus = BUILD_THREADS,
               time_min = 30,
    shell: "cargo build -j {threads} --release --target-dir '{params.rust_dir}' --manifest-path 'implementation/Cargo.toml' --offline 2>&1 | tee '{log.log}'"

localrules: download_bcalm2_gfa_converter
rule download_bcalm2_gfa_converter:
    output: CONVERT_TO_GFA_BINARY,
    conda: "config/conda-download-env.yml"
    params: external_software_scripts_dir = EXTERNAL_SOFTWARE_SCRIPTS_DIR,
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_scripts_dir}'
        cd '{params.external_software_scripts_dir}'

        rm -rf convertToGFA.py
        wget https://raw.githubusercontent.com/GATB/bcalm/v2.2.3/scripts/convertToGFA.py
        chmod u+x convertToGFA.py
        """

localrules: install_contig_validator
rule install_contig_validator:
    input:  sdsl = SDSL_DIR,
    output: dir = directory(CONTIG_VALIDATOR_DIR),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-contigvalidator-env.yml"
    threads: MAX_THREADS
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf ContigValidator
        git clone --recursive https://github.com/mayankpahadia1993/ContigValidator.git
        cd ContigValidator/src
        echo 'count_kmers: count_kmers_kmc' >> Makefile
        sed -i 's\\count_kmers: count_kmers_kmc.cpp KMC/kmc_api/kmc_file.o\\count_kmers_kmc: count_kmers_kmc.cpp KMC/kmc_api/kmc_file.o\\g' Makefile
        LIBRARY_PATH="../../sdsl-lite/lib" CPATH="../../sdsl-lite/include" make -j {threads}
        """

localrules: install_quast
rule install_quast:
    output: script = QUAST_BINARY,
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    threads: 1
    shell: """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf quast
        git clone https://github.com/sebschmi/quast
        cd quast
        git checkout 618b469acef6921da1107022193def891273f6a6
    """

localrules: install_sdsl
rule install_sdsl:
    output: dir = SDSL_DIR,
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-contigvalidator-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf sdsl-lite
        git clone https://github.com/simongog/sdsl-lite.git
        cd sdsl-lite
        git checkout v2.1.1
        HOME=`pwd` ./install.sh
        """

localrules: install_ratatosk
rule install_ratatosk:
    output: binary = RATATOSK_BINARY,
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda: "config/conda-install-ratatosk-env.yml"
    threads: 1
    shell: """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf Ratatosk
        git clone --recursive https://github.com/GuillaumeHolley/Ratatosk.git
        cd Ratatosk
        git checkout --recurse-submodules 74ca617afb20a7c24d73d20f2dcdf223db303496

        mkdir build
        cd build
        cmake ..
        make -j {threads}
        """

localrules: download_wtdbg2
rule download_wtdbg2:
    output: makefile = os.path.join(WTDBG2_DIR, "Makefile"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf wtdbg2
        git clone https://github.com/sebschmi/wtdbg2.git
        cd wtdbg2
        git checkout 78c3077b713aaee48b6c0835105ce6c666f6e796

        sed -i 's:CFLAGS=:CFLAGS=-I${{CONDA_PREFIX}}/include -L${{CONDA_PREFIX}}/lib :g' Makefile
        """

rule build_wtdbg2:
    input:  makefile = os.path.join(WTDBG2_DIR, "Makefile"),
    output: kbm2 = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "wtdbg2/kbm2"),
            pgzf = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "wtdbg2/pgzf"),
            wtdbg2 = WTDBG2_BINARY,
            wtdbg_cns = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "wtdbg2/wtdbg-cns"),
            wtpoa_cns = WTDBG2_CONSENSUS_BINARY,
    params: wtdbg2_dir = WTDBG2_DIR,
    conda: "config/conda-install-wtdbg2-env.yml"
    threads: BUILD_THREADS,
    resources:
        cpus = BUILD_THREADS,
    shell: """
        cd '{params.wtdbg2_dir}'
        make CC=x86_64-conda-linux-gnu-gcc -j {threads}
        """

rule install_sim_it:
    output: SIM_IT_BINARY,
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-download-env.yml"
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        wget -O sim-it.tar.gz https://github.com/ndierckx/Sim-it/archive/refs/tags/Sim-it1.2.tar.gz

        rm -rf Sim-it-Sim-it1.2
        rm -rf sim-it
        tar -xf sim-it.tar.gz
        mv Sim-it-Sim-it1.2/ sim-it/
        mv sim-it/Sim-it1.2.pl sim-it/sim-it.pl
    """

localrules: download_flye
rule download_flye:
    output: flye_marker = os.path.join(FLYE_DIR, ".git", "HEAD"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf Flye
        git clone https://github.com/sebschmi/Flye
        cd Flye
        git checkout 38921327d6c5e57a59e71a7181995f2f0c04be75

        mv bin/flye bin/flye.disabled # rename such that snakemake does not delete it
        """

# Do not make localrule, ensure it is compiled on the correct CPU.
# Otherwise, the compiler might generate unsupported instructions.
rule build_flye:
    input:  flye_marker = os.path.join(FLYE_DIR, ".git", "HEAD"),
    params: flye_directory = FLYE_DIR,
    output: script = FLYE_BINARY,
    conda:  "config/conda-install-flye-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.flye_directory}'

        export CXX=x86_64-conda-linux-gnu-g++
        export CC=x86_64-conda-linux-gnu-gcc
        # export INCLUDES=-I/usr/include/ # Somehow this is not seen by minimap's Makefile, so we had to change it in our custom version of Flye
        # The following also doesn't seem to work when building minimap, so again we had to modify minimap's Makefile
        # export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:=''}} # Redirect library path to include conda libraries
        # make # This does not create the python script anymore

        /usr/bin/env python3 setup.py install

        mv bin/flye.disabled bin/flye # was renamed such that snakemake does not delete it
        """

localrules: download_mdbg
rule download_mdbg:
    output: mdbg_cargo_toml = MDBG_CARGO_TOML,
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
            mdbg_target_directory = os.path.join(MDBG_DIR, "target"),
    conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf rust-mdbg
        git clone https://github.com/sebschmi/rust-mdbg
        cd rust-mdbg
        git checkout 4ff0122a8c63210820ba0341fa7365d6ac216612

        cargo fetch
        
        # rename such that snakemake does not delete them
        mv utils/magic_simplify utils/magic_simplify.original
        mv utils/multik utils/multik.original
        """

rule build_mdbg:
    input:  mdbg_cargo_toml = MDBG_CARGO_TOML,
    output: binary = MDBG_BINARY,
            simplify = MDBG_SIMPLIFY,
            multi_k = MDBG_MULTI_K,
    params: mdbg_directory = MDBG_DIR,
            mdbg_target_directory = os.path.abspath(os.path.join(MDBG_DIR, "target")),
            rust_mdbg = os.path.abspath(os.path.join(MDBG_DIR, "target", "release", "rust-mdbg")),
            to_basespace = os.path.abspath(os.path.join(MDBG_DIR, "target", "release", "to_basespace")),
    conda:  "config/conda-install-mdbg-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.mdbg_directory}'
        cargo --offline build --release -j {threads} --target-dir '{params.mdbg_target_directory}'
        
        # were renamed such that snakemake does not delete them
        cp utils/magic_simplify.original utils/magic_simplify
        cp utils/multik.original utils/multik

        # use built binaries instead of rerunning cargo
        sed -i 's:cargo run --manifest-path .DIR/../Cargo.toml --release:'"'"'{params.rust_mdbg}'"'"':g' utils/multik
        sed -i 's:cargo run --manifest-path .DIR/../Cargo.toml --release --bin to_basespace --:'"'"'{params.to_basespace}'"'"':g' utils/magic_simplify
        """

localrules: download_lja
rule download_lja:
    output: lja_marker = os.path.join(LJA_DIR, "CMakeLists.txt"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-download-env.yml"
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf LJA
        git clone https://github.com/AntonBankevich/LJA
        cd LJA
        git checkout 99f93262c50ff269ee28707f7c3bb77ea00eb576

        #sed -i 's/find_package(OpenMP)//g' CMakeLists.txt
        #sed -i "s:\${{OpenMP_CXX_FLAGS}}:-L${{CONDA_PREFIX}}/lib -lgomp :g" CMakeLists.txt
        #sed -i "s:\${{OpenMP_C_FLAGS}}:-L${{CONDA_PREFIX}}/lib -lgomp :g" CMakeLists.txt
        #sed -i "s:\${{OpenMP_EXE_LINKER_FLAGS}}:-L${{CONDA_PREFIX}}/lib -lgomp :g" CMakeLists.txt
        """

rule build_lja:
    input:  lja_marker = os.path.join(LJA_DIR, "CMakeLists.txt"),
    output: binary = LJA_BINARY,
    params: lja_directory = LJA_DIR,
    conda:  "config/conda-install-lja-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.lja_directory}'

        export CXX=x86_64-conda-linux-gnu-g++
        export CC=x86_64-conda-linux-gnu-gcc

        cmake .
        make -j {threads}
        """

localrules: download_hifiasm
rule download_hifiasm:
    output: hifiasm_marker = os.path.join(HIFIASM_DIR, "Makefile"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-download-env.yml"
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf hifiasm
        git clone https://github.com/sebschmi/hifiasm
        cd hifiasm
        git checkout c914c80547d8cdcfef392291831d6b2fb3b011f5
        """

rule build_hifiasm:
    input:  hifiasm_marker = os.path.join(HIFIASM_DIR, "Makefile"),
    output: binary = HIFIASM_BINARY,
    params: hifiasm_directory = HIFIASM_DIR,
    conda:  "config/conda-install-hifiasm-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.hifiasm_directory}'

        make CXX=x86_64-conda-linux-gnu-g++ CC=x86_64-conda-linux-gnu-gcc CXXFLAGS=-I${{CONDA_PREFIX}}/include -j {threads}
        #make -j {threads}
        """

localrules: download_homopolymer_compress_rs
rule download_homopolymer_compress_rs:
    output: cargo_toml = os.path.join(HOMOPOLYMER_COMPRESS_RS_DIR, "Cargo.toml"),
    log:    log = os.path.join(HOMOPOLYMER_COMPRESS_RS_DIR, "download.log"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf homopolymer-compress-rs
        git clone https://github.com/sebschmi/homopolymer-compress-rs.git
        cd homopolymer-compress-rs
        git checkout 9a979197d2c762f03442a5d584d8c849c9f5ea8e

        cargo fetch
        """

rule build_homopolymer_compress_rs:
    input:  cargo_toml = os.path.join(HOMOPOLYMER_COMPRESS_RS_DIR, "Cargo.toml"),
    output: binary = HOMOPOLYMER_COMPRESS_RS_BINARY
    log:    log = os.path.join(HOMOPOLYMER_COMPRESS_RS_DIR, "build.log"),
    params: homopolymer_compress_rs_dir = HOMOPOLYMER_COMPRESS_RS_DIR,
    conda:  "config/conda-rust-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.homopolymer_compress_rs_dir}'
        cargo build --offline --release -j {threads}
        """

localrules: download_wtdbg2_homopolymer_decompression
rule download_wtdbg2_homopolymer_decompression:
    output: cargo_toml = os.path.join(WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR, "Cargo.toml"),
    log:    log = os.path.join(WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR, "download.log"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf wtdbg2-homopolymer-decompression
        git clone https://github.com/sebschmi/wtdbg2-homopolymer-decompression.git
        cd wtdbg2-homopolymer-decompression
        git checkout 3bec6c0b751a70d53312b359171b9a576f67ebb6

        cargo fetch
        """

rule build_wtdbg2_homopolymer_decompression:
    input:  cargo_toml = os.path.join(WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR, "Cargo.toml"),
    output: binary = WTDBG2_HOMOPOLYMER_DECOMPRESSION_BINARY
    log:    log = os.path.join(WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR, "build.log"),
    params: wtdbg2_homopolymer_decompression_dir = WTDBG2_HOMOPOLYMER_DECOMPRESSION_DIR,
    conda:  "config/conda-rust-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.wtdbg2_homopolymer_decompression_dir}'
        cargo build --offline --release -j {threads}
        """

localrules: download_minimap2_homopolymer_decompression
rule download_minimap2_homopolymer_decompression:
    output: cargo_toml = os.path.join(MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR, "Cargo.toml"),
    log:    log = os.path.join(MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR, "download.log"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf minimap2-homopolymer-decompression
        git clone https://github.com/sebschmi/minimap2-homopolymer-decompression.git
        cd minimap2-homopolymer-decompression
        git checkout d9c0b1e1bdf514d519d25b78f9bb6dda9560e375

        cargo fetch
        """

rule build_minimap2_homopolymer_decompression:
    input:  cargo_toml = os.path.join(MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR, "Cargo.toml"),
    output: binary = MINIMAP2_HOMOPOLYMER_DECOMPRESSION_BINARY
    log:    log = os.path.join(MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR, "build.log"),
    params: minimap2_homopolymer_decompression_dir = MINIMAP2_HOMOPOLYMER_DECOMPRESSION_DIR,
    conda:  "config/conda-rust-env.yml"
    threads: BUILD_THREADS
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.minimap2_homopolymer_decompression_dir}'
        cargo build --offline --release -j {threads}
        """

localrules: download_hisim
rule download_hisim:
    output: makefile = os.path.join(HISIM_DIR, "Makefile"),
    log:    os.path.join(HISIM_DIR, "download.log"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-download-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf HI.SIM
        git clone https://github.com/sebschmi/HI.SIM.git
        cd HI.SIM
        git checkout 734c25c4df3775761ca8920a7d2d57dc44cac09c

        sed -i 's:CFLAGS = :CFLAGS = -I${{CONDA_PREFIX}}/include -L${{CONDA_PREFIX}}/lib :g' Makefile
        """

rule build_hisim:
    input:  makefile = os.path.join(HISIM_DIR, "Makefile"),
    output: binaries = [HISIM_MODEL_BINARY, HISIM_SIM_BINARY, HISIM_FASTA_BINARY],
    log:    os.path.join(HISIM_DIR, "build.log"),
    params: hisim_dir = HISIM_DIR,
    conda:  "config/conda-install-hisim-env.yml"
    threads: BUILD_THREADS,
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.hisim_dir}'
        make CC=x86_64-conda-linux-gnu-gcc -j {threads} all
        """

localrules: download_fastk
rule download_fastk:
    output: makefile = os.path.join(FASTK_DIR, "Makefile"),
    log:    os.path.join(FASTK_DIR, "download.log"),
    params: external_software_dir = EXTERNAL_SOFTWARE_ROOTDIR,
    conda:  "config/conda-download-env.yml"
    threads: 1
    shell:  """
        mkdir -p '{params.external_software_dir}'
        cd '{params.external_software_dir}'

        rm -rf FASTK
        git clone https://github.com/thegenemyers/FASTK.git
        cd FASTK
        git checkout 4604bfcdfd9251d05b27fbd5aef38187e9a9c9ad

        sed -i 's:CFLAGS = :CFLAGS = -I${{CONDA_PREFIX}}/include -L${{CONDA_PREFIX}}/lib :g' Makefile
        sed -i 's:CFLAGS   = :CFLAGS = -I${{CONDA_PREFIX}}/include :g' HTSLIB/Makefile
        sed -i 's:LDFLAGS  = :LDFLAGS = -L${{CONDA_PREFIX}}/lib :g' HTSLIB/Makefile
        """

rule build_fastk:
    input:  makefile = os.path.join(FASTK_DIR, "Makefile"),
    output: binaries = [FASTK_BINARY, FASTK_SYMMEX_BINARY],
    log:    os.path.join(FASTK_DIR, "build.log"),
    params: fastk_dir = FASTK_DIR,
    conda:  "config/conda-install-fastk-env.yml"
    threads: BUILD_THREADS,
    resources:
        cpus = BUILD_THREADS,
    shell:  """
        cd '{params.fastk_dir}'
        make CC=x86_64-conda-linux-gnu-gcc -j {threads} deflate.lib
        make CC=x86_64-conda-linux-gnu-gcc -j {threads} libhts.a
        make CC=x86_64-conda-linux-gnu-gcc -j {threads} all
        """

###################################
###### Download requirements ######
###################################

localrules: compile_all
rule compile_all:
    input:  [WTDBG2_HOMOPOLYMER_DECOMPRESSION_BINARY, HOMOPOLYMER_COMPRESS_RS_BINARY, HIFIASM_BINARY, LJA_BINARY, MDBG_BINARY, FLYE_BINARY, WTDBG2_BINARY, QUAST_BINARY, RUST_BINARY],

localrules: download_and_prepare
rule download_and_prepare:
    input:  reads = expand(GENOME_READS, genome = genomes.keys(), read_source = "real", read_simulation_model_source = "none", homopolymer_compression = "none", read_downsampling_factor = "none", uniquify_ids = "no"),
            # correction_reads = expand(CORRECTION_SHORT_READS_FORMAT, corrected_genome = corrected_genomes.keys()),
            references = expand(GENOME_REFERENCE, genome = genomes.keys(), homopolymer_compression = "none", filter_nw = "no", retain_cm = "no", filter_plasmids = "no"),
            quast = QUAST_BINARY,
            wtdbg2 = WTDBG2_BINARY,
            rust = RUST_BINARY,
            ratatosk = RATATOSK_BINARY,

HUMAN_GENOMES = ["HG002_HiFi_20kb_16x", "HG002_HiFi_15kb_37x", "HG002_HiFi_13.5kb_29x"]
localrules: download_human_data
rule download_human_data:
    input:  reads = [GENOME_READS.format(genome = genome, read_source = "real", read_simulation_model_source = "none", homopolymer_compression = "none", read_downsampling_factor = "none", uniquify_ids = "no") for genome in HUMAN_GENOMES],
            references = [GENOME_REFERENCE.format(genome = genome, homopolymer_compression = "none", filter_nw = "no", retain_cm = "no", filter_plasmids = "no") for genome in HUMAN_GENOMES],

##############################
###### Download results ######
##############################

rule sync_turso_results:
    conda: "config/conda-rsync-env.yml"
    shell: """
        mkdir -p data/reports
        rsync --verbose --recursive --no-relative --include="*/" --include="report.pdf" --include="aggregated-report.pdf" --exclude="*" turso:'/proj/sebschmi/git/practical-omnitigs/data/reports/' data/reports
        """

rule sync_tammi_results:
    conda: "config/conda-rsync-env.yml"
    shell: """
        mkdir -p data/reports
        rsync --verbose --recursive --no-relative --include="*/" --include="report.pdf" --include="aggregated-report.pdf" --exclude="*" tammi:'/abga/work/sebschmi/practical-omnitigs/reports/' data/reports
        """
