

rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "results/logs/get-genome.log",
    params:
        ref_file=config["ref"]["file"],
    cache: True
    shell:
        "( if [[ {params.ref_file} == *.gz ]]; then zcat {params.ref_file} > {output}; else cat {params.ref_file} > {output}; fi ) 1>{log} 2>&1"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "results/logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.23.4/bio/samtools/faidx"


rule bwa_mem2_index:
    input:
        "resources/genome.fasta",
    output:
        multiext(
            "resources/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"
        ),
    log:
        "results/logs/bwa-mem2_index.log",
    cache: True
    wrapper:
        "v1.23.4/bio/bwa-mem2/index"


