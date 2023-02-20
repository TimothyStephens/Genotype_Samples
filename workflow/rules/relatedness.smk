

rule vcftools_relatedness2:
    input:
        "results/merged_calls/all.vcf.gz",
    output:
        "results/relatedness/all.vcf.relatedness2"
    log:
        "results/logs/relatedness/vcftools_relatedness2.log"
    params:
        extra=config["vcftools_relatedness2"]["params"],
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools "
        "--gzvcf {input} "
        "--out results/relatedness/all.vcf "
        "{params.extra} "
        "--relatedness2 "
        "1>{log} 2>&1"


rule vcf_clone_detect:
    input:
        "results/merged_calls/all.vcf.gz",
    output:
        "results/relatedness/all.vcf.allelic_similarity.csv"
    log:
        "results/logs/relatedness/vcftools_relatedness2.log"
    params:
        threshold=config["vcf_clone_detect"]["threshold"],
        extra=config["vcf_clone_detect"]["params"],
    conda:
        "../envs/vcf_clone_detect.yaml"
    shell:
        "python3 workflow/scripts/vcf_clone_detect.py "
        "--vcf {input} "
        "--output {output} "
        "--threshold {params.threshold} "
        "{params.extra} "
        "1>{log} 2>&1"




