configfile: "config/default.yaml"

if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

experiment_file_urls = config['experiment_file_urls']
experiment_files = []
for url in experiment_file_urls:
    experiment_files += ["data/" + url.split("/")[-1].split("?")[0]]

rule selftest:
    conda: "config/conda-selftest-env.yaml"
    shell: "conda --version; wget --version"

rule test:
    input: ["data/" + config['test_file'] + ".unitigs.fa.verify", "data/" + config['test_file'] + ".unitigs.fa.deterministic"]
    conda: "config/conda-rust-env.yaml"
    shell: "cmp --silent {input[0]} {input[1]}; cd implementation; cargo test --target-dir '../data/target'"

rule make_bcalm_output_deterministic:
    input: "data/{file}.unitigs.fa"
    output: "data/{file}.unitigs.fa.deterministic"
    shell: "scripts/make_bcalm_output_deterministic.sh '{input}' '{output}'"

rule verify:
    input: "data/{file}.unitigs.fa"
    output: ["data/{file}.unitigs.fa.verify", "data/{file}.unitigs.fa.properties"]
    conda: "config/conda-rust-env.yaml"
    shell: "cd implementation; cargo run --release --target-dir '../data/target' -- --input '../{input}' --output '../{output[0]}' verify > '../{output[1]}'"

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
    output: experiment_files
    params: experiment_file_urls
    shell: "cd data; wget {params}"