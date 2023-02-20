

rule bwa_mem2_map_reads:
    input:
        reads=get_trimmed_reads,
        idx=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        temp("results/mapped/{sample}-{unit}.sorted.bam"),
        temp("results/mapped/{sample}-{unit}.sorted.bam.csi"),
    log:
        "results/logs/bwa_mem/{sample}-{unit}.log",
    params:
        extra=config["bwa_mem2"]["extra"],
        sort="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=config["bwa_mem2"]["sort_extra"] + " --write-index",  # Extra args for samtools/picard.
    threads: config["bwa_mem2"]["threads"]
    wrapper:
        "v1.23.4/bio/bwa-mem2/mem"


rule samtools_merge:
    input:
        lambda w: expand(
            "results/mapped/{sample}-{unit}.sorted.bam",
            sample=w.sample,
            unit=samples.loc[w.sample].unit,
        ),
    output:
        bam="results/mapped/{sample}.bam",
        idx="results/mapped/{sample}.bam.csi",
    log:
        "results/logs/samtools_merge/{sample}.log",
    params:
        extra=config["samtools_merge"]["params"],  # optional additional parameters, excluding --write-index which is implied by idx
    threads: config["samtools_merge"]["threads"]  # Samtools takes additional threads through its option -@
    wrapper:
        "v1.23.4/bio/samtools/merge"


