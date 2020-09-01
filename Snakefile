configfile: "config/default.yml"

if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

experiment_file_urls = config['experiment_file_urls']
experiment_files = []
experiment_file_url_map = dict()
# file names without endings
experiment_file_names = []

for url in experiment_file_urls:
    file_name_with_ending = url.split("/")[-1].split("?")[0]
    experiment_files += ["data/" + file_name_with_ending]
    experiment_file_url_map["data/" + file_name_with_ending] = url

    file_name_without_ending = file_name_with_ending
    while file_name_without_ending.endswith(('.gz', '.fna')):
        file_name_without_ending = file_name_without_ending[:file_name_without_ending.rfind('.')]

    experiment_file_names += ["data/" + file_name_without_ending]

import pathlib, itertools

rust_sources = list(map(str, itertools.chain(pathlib.Path('implementation').glob('**/Cargo.toml'), pathlib.Path('implementation').glob('**/*.rs'))))


rule selftest:
    conda: "config/conda-selftest-env.yml"
    shell: "echo \"snakemake $(snakemake --version)\"; conda --version; wget --version"

rule test_all:
    input: expand("{file}.is_tested", file = experiment_file_names)

rule test:
    input: verify = "data/" + config['test_file'] + ".is_tested"

rule test_single_file:
    input: verify = "{file}.unitigs.fa.verify",
           deterministic = "{file}.unitigs.fa.deterministic"
    output: touch("{file}.is_tested")
    shell: "cmp --silent {input.verify} {input.deterministic}"

rule make_bcalm_output_deterministic:
    input: file = "data/{file}.unitigs.fa", script = "scripts/make_bcalm_output_deterministic.py"
    output: "data/{file}.unitigs.fa.deterministic"
    shell: "python scripts/make_bcalm_output_deterministic.py '{input.file}' '{output}'"

rule verify_genome_graph:
    input: file = "data/{file}.unitigs.fa", binary = "data/target/release/cli"
    output: ["data/{file}.unitigs.fa.verify", "data/{file}.unitigs.fa.properties"]
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' verify --kmer-size 51 --output '{output[0]}' 2>&1 | tee '{output[1]}.tmp' && mv '{output[1]}.tmp' '{output[1]}'"

rule verify_genome:
    input: file = "{file}.fna", binary = ancient("data/target/release/cli")
    output: log = "{file}.is_genome_verified"
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' verify-genome 2>&1 | tee '{output.log}'"

rule circularise_genome:
    input: file = "{file}.fna", binary = ancient("data/target/release/cli")
    output: data = "{file}.fna-circularised", log = "{file}.fna-circularised-log"
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' circularise-genome 2>&1 --output '{output.data}' | tee '{output.log}'"

rule build_rust_release:
    input: "data/rust.is_tested"
    output: "data/target/release/cli"
    shell: "RUSTFLAGS=\"-C target-cpu=native\" cargo build --release --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

rule test_rust:
    input: expand("{source}", source = list(rust_sources))
    output: touch("data/rust.is_tested")
    shell: "cargo test --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

rule bcalm2:
    input: genome = "data/{file}.fna-circularised", verification = "data/{file}.is_genome_verified"
    output: "data/{file}.unitigs.fa"
    conda: "config/conda-bcalm2-env.yml"
    shell: "cd data; bcalm -in ../{input.genome} -kmer-size 51 -abundance-min 1"

rule extract:
    input: "data/{file}.gz"
    output: "data/{file}"
    conda: "config/conda-extract-env.yml"
    shell: "cd data; gunzip -k {wildcards.file}.gz"

rule download_experiment_file:
    output: "{file}.fna.gz"
    params: url = lambda wildcards, output: experiment_file_url_map[str(output)]
    shell: "cd data; wget {params.url}"

rule latex:
    input: "{file}.tex"
    output: "{file}.pdf"
    conda: "config/conda-latex-env.yml"
    shell: "tectonic {input}"

########################
###### Validation ######
########################

rule create_single_validation_tex:
    input: unitigs_contigvalidator = "{file}.unitigs.contigvalidator", unitigs_quast = directory("{file}.unitigs.quast"), script = "scripts/convert_validation_outputs_to_latex.py"
    output: "{file}.unitigs.tex"
    shell: "scripts/convert_validation_outputs_to_latex.py '{input.unitigs_contigvalidator}' '{input.unitigs_quast}/report.tex' {output}"

rule validate_single_file:
    input: "{file}.unitigs.pdf"
    output: touch("{file}.is_validated")

rule validate_all:
    input: expand("{file}.is_validated", file = experiment_file_names)

rule test_validate:
    input: verify = "data/" + config['test_file'] + ".is_validated"

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

rule run_contig_validator_unitigs:
    input: cv = directory("external-software/ContigValidator"), genome = "{file}.unitigs.fa"
    output: result = "{file}.unitigs.contigvalidator", exact_alignments = temp("{file}.unitigs.fa.exact"),
        bwa_bam = temp("{file}.unitigs.fa.bwa.bam"), bwa_bam_bai = temp("{file}.unitigs.fa.bwa.bam.bai"), fn_kmc_pre = temp("{file}.unitigs.fa.fn.kmc_pre"),
        fn_kmc_suf = temp("{file}.unitigs.fa.fn.kmc_suf"), fp_kmc_pre = temp("{file}.unitigs.fa.fp.kmc_pre"), fp_kmc_suf = temp("{file}.unitigs.fa.fp.kmc_suf"),
        kmc_kmc_pre = temp("{file}.unitigs.fa.kmc.kmc_pre"), kmc_kmc_suf = temp("{file}.unitigs.fa.kmc.kmc_suf"), tp_kmc_pre = temp("{file}.unitigs.fa.tp.kmc_pre"),
        tp_kmc_suf = temp("{file}.unitigs.fa.tp.kmc_suf")
    conda: "config/conda-contigvalidator-env.yml"
    shell:
        """
        cd external-software/ContigValidator
        bash run.sh -suffixsave 0 -abundance-min 1 -kmer-size 51 -r '../../{wildcards.file}.fna-circularised' -a '../../{wildcards.file}.unitigs.contigvalidator' -i '../../{wildcards.file}.unitigs.fa'
        """

###################
###### QUAST ######
###################

rule run_quast_unitigs:
    input: cv = directory("external-software/ContigValidator"), unitigs = "{file}.unitigs.fa", reference = "{file}.fna-circularised"
    output: report = directory("{file}.unitigs.quast")
    conda: "config/conda-quast-env.yml"
    shell: "quast -o {output.report} -r {input.reference} {input.unitigs}"