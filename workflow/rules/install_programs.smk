

rule install_ngsRelate:
    output:
        ".snakemake/install_programs/ngsRelate.done"
    log:
        "results/logs/install_programs/ngsRelate.log",
    conda:
        "../envs/ngsRelate.yaml"
    script:
        "../scripts/install_ngsRelate.sh"


rule install_PCAngsd:
    output:
        ".snakemake/install_programs/PCAngsd.done"
    log:
        "results/logs/install_programs/PCAngsd.log",
    conda:
        "../envs/PCAngsd.yaml"
    script:
        "../scripts/install_PCAngsd.sh"


rule install_nQuire:
    output:
        ".snakemake/install_programs/nQuire.done"
    log:
        "results/logs/install_programs/nQuire.log",
    conda:
        "../envs/nQuire.yaml"
    script:
        "../scripts/install_nQuire.sh"


