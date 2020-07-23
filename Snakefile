configfile: "config/default.yaml"

if 'use_conda' in config and config['use_conda']:
    workflow.use_conda = True

rule selftest:
    conda: "config/conda-selftest-env.yaml"
    shell: "conda --version; wget --version"

rule test:
    input: "data/GCF_000008865.2_ASM886v2_genomic.fna.unitigs.fa"

rule bcalm2:
    input: "data/{file}.fna"
    output: "data/{file}.fna.unitigs.fa"
    conda: "config/conda-bcalm2-env.yaml"
    shell: "cd data; bcalm2 -in {input} -kmer-size 21 -abundance-min 1"

rule download_test_file:
    output: "data/GCF_000008865.2_ASM886v2_genomic.fna.gz"
    conda: "config/conda-download-env.yaml"
    shell: "cd data; wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"

rule extract:
    input: "data/{file}.gz"
    output: "data/{file}"
    conda: "config/conda-extract-env.yaml"
    shell: "cd data; gunzip -k {wildcards.file}.gz"