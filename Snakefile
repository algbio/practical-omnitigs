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
    input: file = "{file}.fna", binary = "data/target/release/cli"
    output: log = "{file}.is_genome_verified"
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' verify-genome 2>&1 | tee '{output.log}'"

rule circularise_genome:
    input: file = "{file}.fna", binary = "data/target/release/cli"
    output: data = "{file}.fna-circularised", log = "{file}.fna-circularised-log"
    conda: "config/conda-rust-env.yml"
    shell: "data/target/release/cli --input '{input.file}' circularise-genome 2>&1 --output '{output.data}' | tee '{output.log}'"

rule build_rust_release:
    input: "data/rust.is_tested"
    output: "data/target/release/cli"
    shell: "cargo build --release --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

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
