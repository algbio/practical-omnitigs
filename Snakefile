import pathlib, itertools

###############################
###### Preprocess Config ######
###############################

configfile: "config/default.yml"

# Allow to configure to use conda from the config file
if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

# Preprocess experiments configuration
experiments = config["experiments"]
tests = config["tests"]

for experiment, config in itertools.chain(experiments.items(), tests.items()):
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

###############################
###### Target Generators ######
###############################

def create_experiment_path(experiment):
    return "data/" + experiment + "/"

def create_report_path(experiment, circularised, k, bcalm2_abundance_min):
    return "data/" + experiment + "/" + ("circular" if circularised else "linear") + ".k" + str(k) + "-a" + str(bcalm2_abundance_min) + "-unitigs.report.pdf"

def _generate_read_sim_targets_(experiment, config):
    path = create_experiment_path(experiment)
    for circularised in config["circularised"]:
        yield path + ("circular" if circularised else "linear")

def _generate_bcalm2_targets_(experiment, config):
    for target in _generate_read_sim_targets_(experiment, config):
        for k in config["k"]:
            for bcalm2_abundance_min in config["bcalm2_abundance_min"]:
                    yield target + ".k" + str(k) + "-a" + str(bcalm2_abundance_min) + "-unitigs"

def _generate_report_targets_(experiment, config):
        for target in _generate_bcalm2_targets_(experiment, config):
            yield target + ".report.pdf"

def _generate_test_targets_(experiment, config):
        for target in _generate_bcalm2_targets_(experiment, config):
            yield target + ".is_tested"

def generate_report_targets():
    for experiment, config in experiments.items():
        for target in _generate_report_targets_(experiment, config):
            yield target

def generate_test_report_targets():
    for experiment, config in tests.items():
        for target in _generate_report_targets_(experiment, config):
            yield target

def generate_test_targets():
    for experiment, config in tests.items():
        for target in _generate_test_targets_(experiment, config):
            yield target

######################################
###### Input Genome Preparation ######
######################################

rule separate_linear_and_circular:
    input: filtered = "data/{dir}/filtered.fna", verified = "data/{dir}/is_genome_verified.log", binary = "data/target/release/cli"
    output: circular = "data/{dir}/circular.fna", linear = "data/{dir}/linear.fna", log = "data/{dir}/separate_linear_and_circular.log"
    conda: "config/conda-rust-env.yml"
    shell: "cp '{input.filtered}' '{output.linear}'; data/target/release/cli --input '{input.filtered}' circularise-genome 2>&1 --output '{output.circular}' | tee '{output.log}'"

rule verify_genome:
    input: file = "data/{dir}/filtered.fna", binary = "data/target/release/cli"
    output: log = "data/{dir}/is_genome_verified.log"
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' verify-genome 2>&1 | tee '{output.log}'"

rule filter_genome:
    input: file = "data/{dir}/raw.fna", binary = "data/target/release/cli"
    output: file = "data/{dir}/filtered.fna", log = "data/{dir}/filtered.log"
    params: retain = lambda wildcards: "--retain '" + experiments[wildcards.dir]["filter_retain"] + "'" if "filter_retain" in experiments[wildcards.dir] else ""
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' filter --output '{output.file}' {params.retain} 2>&1 | tee '{output.log}'"

rule extract:
    input: "data/{dir}/raw.fna.gz"
    output: "data/{dir}/raw.fna"
    conda: "config/conda-extract-env.yml"
    shell: "cd 'data/{wildcards.dir}'; gunzip -k raw.fna.gz"

rule download_experiment_file:
    output: "data/{dir}/raw.fna.gz"
    params: url = lambda wildcards, output: experiments[wildcards.dir]["url"]
    conda: "config/conda-download-env.yml"
    shell: "mkdir -p 'data/{wildcards.dir}'; cd 'data/{wildcards.dir}'; wget -O raw.fna.gz {params.url}"

##################
###### Rust ######
##################

rule build_rust_release:
    input: "data/is_rust_tested.log"
    output: "data/target/release/cli"
    shell: "RUSTFLAGS=\"-C target-cpu=native\" cargo build --release --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

rule test_rust:
    input: expand("{source}", source = list(rust_sources))
    output: touch("data/is_rust_tested.log")
    shell: "cargo test --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml' 2>&1 | tee '{output}'"

#####################
###### Testing ######
#####################

rule test:
    input: generate_test_targets()

rule test_single_file:
    input: verify = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.verify",
           deterministic = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.deterministic"
    output: touch("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.is_tested")
    shell: "cmp --silent {input.verify} {input.deterministic}"

rule make_bcalm_output_deterministic:
    input: file = "data/{dir}/{file}-unitigs.fa", script = "scripts/make_bcalm_output_deterministic.py"
    output: file = "data/{dir}/{file}-unitigs.fa.deterministic"
    shell: "python scripts/make_bcalm_output_deterministic.py '{input.file}' '{output.file}'"

rule verify_genome_graph:
    input: file = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa", binary = "data/target/release/cli"
    output: verification_copy = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.verify", log =  "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.properties", latex = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.graphstatistics"
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' verify --kmer-size {wildcards.k} --output '{output.verification_copy}' --latex '{output.latex}' 2>&1 | tee '{output.log}.tmp' && mv '{output.log}.tmp' '{output.log}'"

rule selftest:
    conda: "config/conda-selftest-env.yml"
    shell: "echo \"snakemake $(snakemake --version)\"; conda --version; wget --version"

#######################################
###### Genome Graph Construction ######
#######################################

rule bcalm2:
    input: genome = "data/{dir}/{file}.fna"
    output: unitigs = "data/{dir}/{file,(circular|linear)}.k{k,[0-9]+}-a{abundance_min,[0-9]+}-unitigs.fa",
    #params: tmp = "data/{dir}/{file,(circular|linear)}.k{k,[0-9]+}-a{abundance_min,[0-9]+}-unitigs.bcalm2-tmp/"
    conda: "config/conda-bcalm2-env.yml"
    shell: 
        """
        bcalm -in '{input.genome}' -out '{output.unitigs}' -kmer-size {wildcards.k} -abundance-min {wildcards.abundance_min}
        mv '{output.unitigs}.unitigs.fa' '{output.unitigs}'
        rm data/{wildcards.dir}/{wildcards.file}.k{wildcards.k}-a{wildcards.abundance_min}-unitigs.*.glue.*
        """

###############################
###### Report Generation ######
###############################

rule latex:
    input: "data/{dir}/{file}.report.tex"
    output: "data/{dir}/{file}.report.pdf"
    conda: "config/conda-latex-env.yml"
    shell: "tectonic {input}"

rule create_single_report_tex:
    input: unitigs_contigvalidator = "data/{dir}/{file}.contigvalidator",
           unitigs_quast = directory("data/{dir}/{file}.quast"),
           unitigs_graphstatistics = "data/{dir}/{file}.graphstatistics",
           untitigs_bandage = "data/{dir}/{file}.bandage.png",
           script = "scripts/convert_validation_outputs_to_latex.py",
    output: "data/{dir}/{file}.report.tex"
    shell: "scripts/convert_validation_outputs_to_latex.py '{input.unitigs_contigvalidator}' '{input.unitigs_quast}/report.tex' '{input.unitigs_graphstatistics}' '{wildcards.file}.bandage.png' '{output}'"

rule report_all:
    input: generate_report_targets()

rule test_report:
    input: generate_test_report_targets()

rule png_to_pdf:
    input: "{file}.png"
    output: "{file}.image.pdf"
    conda: "config/conda-imagemagick-env.yml"
    shell: "convert {input} {output}"

#############################
###### ContigValidator ######
#############################

rule install_sdsl:
    output: dir = directory("external-software/sdsl-lite")
    conda: "config/conda-contigvalidator-env.yml"
    shell:
        """
        cd external-software
        git clone https://github.com/simongog/sdsl-lite.git
        cd sdsl-lite
        git checkout v2.1.1
        HOME=`pwd` ./install.sh
        """

rule install_contig_validator:
    input: sdsl = directory("external-software/sdsl-lite")
    output: dir = directory("external-software/ContigValidator")
    conda: "config/conda-contigvalidator-env.yml"
    shell: 
        """
        cd external-software
        git clone --recursive https://github.com/mayankpahadia1993/ContigValidator.git
        cd ContigValidator/src
        LIBRARY_PATH="../../sdsl-lite/lib" CPATH="../../sdsl-lite/include" make
        """

rule run_contig_validator:
    input: cv = directory("external-software/ContigValidator"),
        reads = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa",
        reference = "data/{dir}/{file}.fna"
    output: result = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.contigvalidator",
        exact_alignments = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.exact"),
        bwa_bam = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.bwa.bam"),
        bwa_bam_bai = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.bwa.bam.bai"),
        fn_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.fn.kmc_pre"),
        fn_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.fn.kmc_suf"),
        fp_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.fp.kmc_pre"),
        fp_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.fp.kmc_suf"),
        kmc_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.kmc.kmc_pre"),
        kmc_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.kmc.kmc_suf"),
        tp_kmc_pre = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.tp.kmc_pre"),
        tp_kmc_suf = temp("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa.tp.kmc_suf")
    conda: "config/conda-contigvalidator-env.yml"
    shell:
        """
        cd external-software/ContigValidator
        bash run.sh -suffixsave 0 -abundance-min 1 -kmer-size {wildcards.k} -r '../../{input.reference}' -a '../../{output.result}' -i '../../{input.reads}'
        """

###################
###### QUAST ######
###################

rule run_quast:
    input: reads = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa",
        reference = "data/{dir}/{file}.fna"
    output: report = directory("data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.quast")
    conda: "config/conda-quast-env.yml"
    shell: "quast -o {output.report} -r {input.reference} {input.reads}"

#####################
###### Bandage ######
#####################

rule download_bcalm2_gfa_converter:
    output: "external-software/scripts/convertToGFA.py"
    conda: "config/conda-download-env.yml"
    shell:
        """
        mkdir -p external-software/scripts
        cd external-software/scripts
        wget https://raw.githubusercontent.com/GATB/bcalm/v2.2.3/scripts/convertToGFA.py
        chmod u+x convertToGFA.py
        """

rule convert_bcalm2_output_to_gfa:
    input: fa = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.fa",
        converter = "external-software/scripts/convertToGFA.py"
    output: gfa = "data/{dir}/{file}.k{k}-a{abundance_min}-unitigs.gfa"
    shell: "external-software/scripts/convertToGFA.py {input.fa} {output.gfa} {wildcards.k}"

rule bandage:
    input: "{file}.gfa"
    output: "{file}.bandage.png"
    conda: "config/conda-bandage-env.yml"
    shell: "Bandage image {input} {output} --width 1000 --height 1000"