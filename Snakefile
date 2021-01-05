import pathlib, itertools, sys

###############################
###### Preprocess Config ######
###############################

configfile: "config/default.yml"

# Allow to configure to use conda from the config file
if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

workflow.global_resources["contigvalidator"] = 1
workflow.global_resources["concorde"] = 1

# Preprocess experiments configuration
experiments_bcalm2 = config["experiments"]["bcalm2"]
experiments_wtdbg2 = config["experiments"]["wtdbg2"]

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

for experiment, config in experiments_wtdbg2.items():
    if not "urls" in config:
        print("Missing url in wtdbg2 experiment")
        sys.exit(1)

# Collect all rust sources
rust_sources = list(map(str, itertools.chain(pathlib.Path('implementation').glob('**/Cargo.toml'), pathlib.Path('implementation').glob('**/*.rs'))))

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
    input: filtered = "data/{dir}/filtered.fna", verified = "data/{dir}/is_genome_verified.log", binary = "data/target/release/cli"
    output: circular = "data/{dir}/circular.fna", linear = "data/{dir}/linear.fna", log = "data/{dir}/separate_linear_and_circular.log"
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
    shell: "mkdir -p 'data/{wildcards.dir}'; cd 'data/{wildcards.dir}'; wget -O raw.fna.gz {params.url}"

######################################
###### wtdbg2 Input Preparation ######
######################################

rule download_wtdbg2_input:
    output: reads = "data/{dir}/reads.fa", reference = "data/{dir}/reference.fa.gz"
    params: url = lambda wildcards, output: "'" + "' '".join(experiments_wtdbg2[wildcards.dir]["urls"]) + "'",
            reference = lambda wildcards, output: experiments_wtdbg2[wildcards.dir]["reference"]
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: "mkdir -p 'data/{wildcards.dir}'; wget -O '{output.reads}' {params.url}; wget -O '{output.reference}' '{params.reference}'"

rule uniquify_fasta_ids:
    input: reads = "data/{dir}/reads.fa", script = "scripts/uniquify_fasta_ids.py"
    output: reads = "data/{dir}/reads.uniqified.fa", log = "data/{dir}/uniquify.log"
    conda: "config/conda-uniquify-env.yml"
    threads: 1
    shell: "python3 '{input.script}' '{input.reads}' '{output.reads}' 2>&1 | tee '{output.log}'"

##################
###### Rust ######
##################

rule build_rust_release:
    input: "data/is_rust_tested.log"
    output: "data/target/release/cli"
    conda: "config/conda-rust-env.yml"
    threads: workflow.cores
    shell: "RUSTFLAGS=\"-C target-cpu=native\" cargo build -j {threads} --release --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

rule test_rust:
    input: expand("{source}", source = list(rust_sources))
    output: touch("data/is_rust_tested.log")
    conda: "config/conda-rust-env.yml"
    threads: workflow.cores
    shell: "cargo test -j {threads} --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml' 2>&1 | tee '{output}'"

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

### wtdbg2 ###

# rule compute_omnitigs_wtdbg2:
#     input: nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot", raw_reads = "data/{dir}/reads.uniqified.fa", binary = "data/target/release/cli"
#     output: file = "data/{dir}/wtdbg2.omnitigs.contigwalks", log = "data/{dir}/wtdbg2.omnitigs.log", latex = "data/{dir}/wtdbg2.omnitigs.tex"
#     threads: 1
#     shell: "'{input.binary}' compute-omnitigs --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

# rule compute_trivial_omnitigs_wtdbg2:
#     input: nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot", raw_reads = "data/{dir}/reads.uniqified.fa", binary = "data/target/release/cli"
#     output: file = "data/{dir}/wtdbg2.trivialomnitigs.contigwalks", log = "data/{dir}/wtdbg2.trivialomnitigs.log", latex = "data/{dir}/wtdbg2.trivialomnitigs.tex"
#     threads: 1
#     shell: "'{input.binary}' compute-trivial-omnitigs --output-as-wtdbg2-node-ids --non-scc --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

# rule compute_unitigs_wtdbg2:
#     input: nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot", ctg_lay = "data/{dir}/wtdbg2.wtdbg2.ctg.lay", raw_reads = "data/{dir}/reads.uniqified.fa", binary = "data/target/release/cli"
#     output: file = "data/{dir}/wtdbg2.unitigs.contigwalks", log = "data/{dir}/wtdbg2.unitigs.log", latex = "data/{dir}/wtdbg2.unitigs.tex"
#     threads: 1
#     shell: "'{input.binary}' compute-unitigs --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

### injected wtdbg2 ###

rule compute_injectable_omnitigs_wtdbg2:
    input: nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot", raw_reads = "data/{dir}/reads.uniqified.fa", binary = "data/target/release/cli"
    output: file = "data/{dir}/wtdbg2.injected-omnitigs{algo_suffix}.contigwalks", log = "data/{dir}/wtdbg2.injected-omnitigs{algo_suffix}.log", latex = "data/{dir}/wtdbg2.injected-omnitigs{algo_suffix}.tex"
    wildcard_constraints: algo_suffix = ".*"
    threads: 1
    shell: "'{input.binary}' compute-omnitigs --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

rule compute_injectable_trivial_omnitigs_wtdbg2:
    input: nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot", raw_reads = "data/{dir}/reads.uniqified.fa", binary = "data/target/release/cli"
    output: file = "data/{dir}/wtdbg2.injected-trivialomnitigs{algo_suffix}.contigwalks", log = "data/{dir}/wtdbg2.injected-trivialomnitigs{algo_suffix}.log", latex = "data/{dir}/wtdbg2.injected-trivialomnitigs{algo_suffix}.tex"
    wildcard_constraints: algo_suffix = ".*"
    threads: 1
    shell: "'{input.binary}' compute-trivial-omnitigs --output-as-wtdbg2-node-ids --non-scc --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

rule compute_injectable_unitigs_wtdbg2:
    input: nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot", ctg_lay = "data/{dir}/wtdbg2.wtdbg2.ctg.lay", raw_reads = "data/{dir}/reads.uniqified.fa", binary = "data/target/release/cli"
    output: file = "data/{dir}/wtdbg2.injected-unitigs{algo_suffix}.contigwalks", log = "data/{dir}/wtdbg2.injected-unitigs{algo_suffix}.log", latex = "data/{dir}/wtdbg2.injected-unitigs{algo_suffix}.tex"
    wildcard_constraints: algo_suffix = ".*"
    threads: 1
    shell: "'{input.binary}' compute-unitigs --output-as-wtdbg2-node-ids --file-format wtdbg2 --input '{input.nodes}' --input '{input.reads}' --input '{input.raw_reads}' --input '{input.dot}' --output '{output.file}' --latex '{output.latex}' 2>&1 | tee '{output.log}'"

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

####################
###### wtdbg2 ######
####################

rule install_wtdbg2:
    output: kbm2 = "external-software/wtdbg2/kbm2", pgzf = "external-software/wtdbg2/pgzf", wtdbg2 = "external-software/wtdbg2/wtdbg2", wtdbg_cns = "external-software/wtdbg2/wtdbg-cns", wtpoa_cns = "external-software/wtdbg2/wtpoa-cns"
    conda: "config/conda-download-env.yml"
    threads: 1
    shell: """
    mkdir -p external-software
    cd external-software

    git clone https://github.com/sebschmi/wtdbg2.git
    cd wtdbg2
    git checkout 11011ee9c27f10d08fd302bbb68f4a6aba5f9748
    make
    """

rule wtdbg2_complete:
    input: reads = "data/{dir}/reads.uniqified.fa", binary = "external-software/wtdbg2/wtdbg2"
    output: original_nodes = "data/{dir}/wtdbg2.wtdbg2.1.nodes", nodes = "data/{dir}/wtdbg2.wtdbg2.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2.3.dot.gz", clips = "data/{dir}/wtdbg2.wtdbg2.clps", kbm = "data/{dir}/wtdbg2.wtdbg2.kbm", ctg_lay = "data/{dir}/wtdbg2.wtdbg2.ctg.lay.gz", log = "data/{dir}/wtdbg2.wtdbg2.log"
    threads: workflow.cores
    shell: "{input.binary} -x rs -g 100m -i '{input.reads}' -t {threads} -fo 'data/{wildcards.dir}/wtdbg2.wtdbg2' --dump-kbm '{output.kbm}' 2>&1 | tee '{output.log}'"

rule wtdbg2_complete_without_fragment_assembly:
    input: reads = "data/{dir}/reads.uniqified.fa", binary = "external-software/wtdbg2/wtdbg2"
    output: original_nodes = "data/{dir}/wtdbg2.wtdbg2-sfa.1.nodes", nodes = "data/{dir}/wtdbg2.wtdbg2-sfa.3.nodes", reads = "data/{dir}/wtdbg2.wtdbg2-sfa.3.reads", dot = "data/{dir}/wtdbg2.wtdbg2-sfa.3.dot.gz", clips = "data/{dir}/wtdbg2.wtdbg2-sfa.clps", kbm = "data/{dir}/wtdbg2.wtdbg2-sfa.kbm", ctg_lay = "data/{dir}/wtdbg2.wtdbg2-sfa.ctg.lay.gz", log = "data/{dir}/wtdbg2.wtdbg2-sfa.log"
    threads: workflow.cores
    shell: "{input.binary} -x rs -g 100m --skip-fragment-assembly -i '{input.reads}' -t {threads} -fo 'data/{wildcards.dir}/wtdbg2.wtdbg2-sfa' --dump-kbm '{output.kbm}' 2>&1 | tee '{output.log}'"

rule wtdbg2_inject_contigs:
    input: reads = "data/{dir}/reads.uniqified.fa", contigs = "data/{dir}/wtdbg2.injected-{algorithm}.contigwalks", clips = "data/{dir}/wtdbg2.wtdbg2.clps", nodes = "data/{dir}/wtdbg2.wtdbg2.1.nodes", kbm = "data/{dir}/wtdbg2.wtdbg2.kbm", binary = "external-software/wtdbg2/wtdbg2"
    output: ctg_lay = "data/{dir}/wtdbg2.injected-{algorithm}.ctg.lay.gz", log = "data/{dir}/wtdbg2.injected-{algorithm}.log"
    params: genome_length = lambda wildcards, output: "-g 100m", #if wildcards.algorithm == "unitigs" else ""
            skip_fragment_assembly = lambda wildcards, output: "--skip-fragment-assembly" if "-sfa" in wildcards.algorithm else ""
    threads: workflow.cores
    shell: "{input.binary} -x rs {params.genome_length} {params.skip_fragment_assembly} -i '{input.reads}' -t {threads} -fo 'data/{wildcards.dir}/wtdbg2.injected-{wildcards.algorithm}' --inject-unitigs '{input.contigs}' --load-nodes '{input.nodes}' --load-clips '{input.clips}' --load-kbm '{input.kbm}' 2>&1 | tee '{output.log}'"

rule wtdbg2_consensus:
    input: reads = "data/{dir}/reads.uniqified.fa", contigs = "data/{dir}/wtdbg2.{algorithm}.ctg.lay", binary = "external-software/wtdbg2/wtpoa-cns"
    output: consensus = "data/{dir}/wtdbg2.{algorithm}.raw.fa"
    threads: workflow.cores
    shell: "{input.binary} -t {threads} -i '{input.contigs}' -fo '{output.consensus}'"

###############################
###### Report Generation ######
###############################

rule latex:
    input: "data/{dir}/{file}report.tex"
    output: "data/{dir}/{file}report.pdf"
    conda: "config/conda-latex-env.yml"
    threads: 1
    shell: "tectonic {input}"

rule create_single_report_tex:
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
    threads: 1
    shell: "scripts/convert_validation_outputs_to_latex.py '{input.genome_name}' '{input.graphstatistics}' '{input.bcalm2_bandage}' 'none' '{output}' uni '{params.prefix}.unitigs' 'Y-to-V' '{params.prefix}.trivialomnitigs' omni '{params.prefix}.omnitigs'"

rule create_single_wtdbg2_report_tex:
    input: unitigs = "data/{dir}/wtdbg2.injected-unitigs-sfa.tex",
           trivialomnitigs = "data/{dir}/wtdbg2.injected-trivialomnitigs-sfa.tex",
           #omnitigs_quast = "data/{dir}/wtdbg2.injected-omnitigs-sfa.tex",
           injected_unitigs = "data/{dir}/wtdbg2.injected-unitigs.tex",
           injected_trivialomnitigs = "data/{dir}/wtdbg2.injected-trivialomnitigs.tex",
           #injected_omnitigs_quast = "data/{dir}/wtdbg2.injected-omnitigs.tex",
           unitigs_quast = "data/{dir}/wtdbg2.injected-unitigs-sfa.quast",
           trivialomnitigs_quast = "data/{dir}/wtdbg2.injected-trivialomnitigs-sfa.quast",
           #omnitigs_quast = "data/{dir}/wtdbg2.omnitigs.quast",
           injected_unitigs_quast = "data/{dir}/wtdbg2.injected-unitigs.quast",
           injected_trivialomnitigs_quast = "data/{dir}/wtdbg2.injected-trivialomnitigs.quast",
           #injected_omnitigs_quast = "data/{dir}/wtdbg2.injected-omnitigs.quast",
           wtdbg2_quast = "data/{dir}/wtdbg2.wtdbg2.quast",
           wtdbg2_sfa_quast = "data/{dir}/wtdbg2.wtdbg2-sfa.quast",
           combined_eaxmax_plot = "data/{dir}/wtdbg2.wtdbg2-eaxmax-plot.pdf",
           script = "scripts/convert_validation_outputs_to_latex.py",
    output: "data/{dir}/wtdbg2.wtdbg2-report.tex"
    params: prefix = "data/{dir}/wtdbg2"
    threads: 1
    shell: """echo '{wildcards.dir}' > 'data/{wildcards.dir}/name.txt'
              '{input.script}' 'data/{wildcards.dir}/name.txt' 'none' 'none' '{input.combined_eaxmax_plot}' '{output}' 'inj uni sfa' '{params.prefix}.injected-unitigs-sfa' 'inj Y-to-V sfa' '{params.prefix}.injected-trivialomnitigs-sfa' 'wtdbg2 sfa' '{params.prefix}.wtdbg2-sfa' 'inj uni' '{params.prefix}.injected-unitigs' 'inj Y-to-V' '{params.prefix}.injected-trivialomnitigs' wtdbg2 '{params.prefix}.wtdbg2'"""
              #scripts/convert_validation_outputs_to_latex.py 'data/{wildcards.dir}/name.txt' 'none' 'none' '{output}' uni '{params.prefix}.unitigs' Y-to-V '{params.prefix}.trivialomnitigs' omni '{params.prefix}.omnitigs' 'inj uni' '{params.prefix}.unitigs' 'inj Y-to-V' '{params.prefix}.trivialomnitigs' 'inj omni' '{params.prefix}.omnitigs' wtdbg2 '{params.prefix}.wtdbg2'"""

rule create_combined_eaxmax_graph:
    input: unitigs_quast = "data/{dir}/wtdbg2.injected-unitigs-sfa.quast",
           trivialomnitigs_quast = "data/{dir}/wtdbg2.injected-trivialomnitigs-sfa.quast",
           #omnitigs_quast = "data/{dir}/wtdbg2.omnitigs.quast",
           injected_unitigs_quast = "data/{dir}/wtdbg2.injected-unitigs.quast",
           injected_trivialomnitigs_quast = "data/{dir}/wtdbg2.injected-trivialomnitigs.quast",
           #injected_omnitigs_quast = "data/{dir}/wtdbg2.injected-omnitigs.quast",
           wtdbg2_quast = "data/{dir}/wtdbg2.wtdbg2.quast",
           wtdbg2_sfa_quast = "data/{dir}/wtdbg2.wtdbg2-sfa.quast",
           script = "scripts/create_combined_eaxmax_plot.py",
    output: "data/{dir}/wtdbg2.wtdbg2-eaxmax-plot.pdf",
    conda: "config/conda-seaborn-env.yml"
    threads: 1
    shell: "python3 '{input.script}' 'data/{wildcards.dir}/' '{output}'"

rule report_all:
    input: generate_report_targets()

rule report_wtdbg2:
    input: generate_wtdbg2_report_targets()

rule test_report:
    input: generate_test_report_targets()

rule png_to_pdf:
    input: "{file}.png"
    output: "{file}.image.pdf"
    conda: "config/conda-imagemagick-env.yml"
    threads: 1
    shell: "convert {input} {output}"

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

###################
###### QUAST ######
###################

rule install_quast:
    output: script = "external-software/quast/quast.py", script_directory = directory("external-software/quast/")
    conda: "config/conda-install-env.yml"
    threads: 1
    shell: """
    mkdir -p external-software
    cd external-software

    git clone https://github.com/sebschmi/quast
    cd quast
    git checkout ef6e48274c40b2603dac63bc31450e70cfc4c690
    """

rule run_quast:
    input: reads = "data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.fa",
        reference = "data/{dir}/{file}.fna",
        script = "external-software/quast/quast.py",
        script_directory = "external-software/quast/"
    output: report = directory("data/{dir}/{file}.k{k}-a{abundance_min}.{algorithm}.quast")
    conda: "config/conda-quast-env.yml"
    threads: workflow.cores / 8
    shell: "{input.script} -t {threads} -o {output.report} -r {input.reference} {input.reads}"

rule run_quast_wtdbg2:
    input: contigs = "data/{dir}/wtdbg2.{algorithm}.raw.fa",
        reference = "data/{dir}/reference.fa",
        script = "external-software/quast/quast.py",
        script_directory = "external-software/quast/"
    output: report = directory("data/{dir}/wtdbg2.{algorithm}.quast")
    conda: "config/conda-quast-env.yml"
    threads: workflow.cores / 8
    shell: "{input.script} -t {threads} --fragmented -o {output.report} -r {input.reference} {input.contigs}"

rule test_quast_wtdbg2:
    input: contigs = "data/C.elegans/wtdbg2.wtdbg2.raw.fa",
        reference = "data/C.elegans/reference.fa",
        script = "../quast/quast.py"
    params: report = "data/C.elegans/wtdbg2.wtdbg2.quast"
    conda: "config/conda-quast-env.yml"
    threads: workflow.cores / 8
    shell: """{input.script} -t {threads} --fragmented -o {params.report} -r {input.reference} {input.contigs}"""

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
