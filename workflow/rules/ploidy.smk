

rule nQuire_create:
    input:
        bam=rules.samtools_merge.output.bam,
        idx=rules.samtools_merge.output.idx,
        programs=rules.install_nQuire.output,
    output:
        "results/ploidy/{sample}.bin",
    log:
        "results/logs/ploidy/nQuire_create.{sample}.log",
    params:
        extra=config["nQuire"]["create_params"],
        min_quality=config["nQuire"]["min_quality"],
        min_coverage=config["nQuire"]["min_coverage"],
    conda:
        "../envs/nQuire.yaml"
    shell:
        "(nQuire create -q {params.min_quality} -c {params.min_coverage} -x -b {input.bam} -o $(echo {output} | sed -e 's/.bin//') {params.extra}) 1>{log} 2>&1"


rule nQuire_denoise:
    input:
        nQbin=rules.nQuire_create.output,
        programs=rules.install_nQuire.output,
    output:
        "results/ploidy/{sample}.denoised.bin",
    log:
        "results/logs/ploidy/nQuire_denoise.{sample}.log",
    params:
        extra=config["nQuire"]["denoise_params"],
    conda:
        "../envs/nQuire.yaml"
    shell:
        "(nQuire denoise -i {input.nQbin} -o $(echo {output} | sed -e 's/.bin//') {params.extra}) 1>{log} 2>&1"


rule nQuire_lrdmodel:
    input:
        nQbins=[
              *expand("results/ploidy/{sample}.bin", sample=samples.sample_id), 
              *expand("results/ploidy/{sample}.denoised.bin", sample=samples.sample_id),
        ],
        programs=rules.install_nQuire.output,
    output:
        "results/ploidy/nQuire_lrdmodel.txt",
    log:
        "results/logs/ploidy/nQuire_lrdmodel.log",
    params:
        extra=config["nQuire"]["lrdmodel_params"],
    threads: config["nQuire"]["threads"]
    conda:
        "../envs/nQuire.yaml"
    shell:
        "(nQuire lrdmodel {params.extra} -t {threads} {input.nQbins} > {output}) 1>{log} 2>&1"


rule format_nQuire_results:
    input:
        "results/ploidy/nQuire_lrdmodel.txt",
    output:
        "results/final/nQuire.tsv",
    log:
        "results/logs/ploidy/format_nQuire_results.log",
    script:
        "../scripts/format_nQuire_results.sh"


