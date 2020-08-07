configfile: "config/default.yaml"

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
    conda: "config/conda-selftest-env.yaml"
    shell: "conda --version; wget --version"

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
    input: "data/{file}.unitigs.fa"
    output: "data/{file}.unitigs.fa.deterministic"
    shell: "scripts/make_bcalm_output_deterministic.sh '{input}' '{output}'"

rule verify:
    input: file = "data/{file}.unitigs.fa", binary = "data/target/release/cli"
    output: ["data/{file}.unitigs.fa.verify", "data/{file}.unitigs.fa.properties"]
    conda: "config/conda-rust-env.yaml"
    shell: "data/target/release/cli --input '{input.file}' --output '{output[0]}' verify 2>&1 | tee '{output[1]}.tmp' && mv '{output[1]}.tmp' '{output[1]}'"

rule build_rust:
    input: expand("{source}", source = list(rust_sources))
    output: "data/target/release/cli"
    shell: "cargo test --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'; cargo build --release --target-dir 'data/target' --manifest-path 'implementation/Cargo.toml'"

rule bcalm2:
    input: "data/{file}.fna"
    output: "data/{file}.unitigs.fa"
    conda: "config/conda-bcalm2-env.yaml"
    shell: "cd data; bcalm -in ../{input} -kmer-size 21 -abundance-min 1"

rule extract:
    input: "data/{file}.gz"
    output: "data/{file}"
    conda: "config/conda-extract-env.yaml"
    shell: "cd data; gunzip -k {wildcards.file}.gz"

rule download_experiment_file:
    output: "{file}.fna.gz"
    params: url = lambda wildcards, output: experiment_file_url_map[str(output)]
    shell: "cd data; wget {params.url}"