

rule deepvariant:
    input:
        bam=rules.samtools_merge.output.bam,
        idx=rules.samtools_merge.output.idx,
        ref=rules.get_genome.output,
        ref_idx=rules.genome_faidx.output,
    output:
        vcf="results/calls/{sample}.vcf.gz",
        report=report(
            "results/calls/{sample}.visual_report.html",
            caption="../report/vcf.rst",
            category="Calls",
        ),
    params:
        model=config["deepvariant"]["model"],
        extra=config["deepvariant"]["extra"],
        tmp_dir="results/calls/{sample}.tmp",
    threads: config["deepvariant"]["threads"]
    log:
        "results/logs/deepvariant/{sample}/stdout.log",
    container:
        "docker://google/deepvariant:1.4.0"
    shell:
        "run_deepvariant "
        "--model_type {params.model} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "--sample_name {wildcards.sample} "
        "--logging_dir results/logs/deepvariant/{wildcards.sample} "
        "--output_vcf {output.vcf} "
        "--intermediate_results_dir {params.tmp_dir} "
        "--num_shards {threads} "
        "1>{log} 2>&1"


rule deepvariant_gvcf:
    input:
        bam=rules.samtools_merge.output.bam,
        idx=rules.samtools_merge.output.idx,
        ref=rules.get_genome.output,
        ref_idx=rules.genome_faidx.output,
    output:
        vcf="results/individual_calls/{sample}.vcf.gz",
        report=report(
            "results/individual_calls/{sample}.visual_report.html",
            caption="../report/vcf.rst",
            category="Calls",
        ),
        gvcf="results/individual_calls/{sample}.g.vcf.gz",
    params:
        model=config["deepvariant_gvcf"]["model"],
        extra=config["deepvariant_gvcf"]["extra"],
        tmp_dir="results/individual_calls/{sample}.tmp",
    threads: config["deepvariant_gvcf"]["threads"]
    log:
        "results/logs/deepvariant_gvcf/{sample}/stdout.log",
    container:
        "docker://google/deepvariant:1.4.0"
    shell:
        "run_deepvariant "
        "--model_type {params.model} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "--sample_name {wildcards.sample} "
        "--logging_dir results/logs/deepvariant/{wildcards.sample} "
        "--output_vcf {output.vcf} "
        "--output_gvcf {output.gvcf} "
        "--intermediate_results_dir {params.tmp_dir} "
        "--num_shards {threads} "
        "1>{log} 2>&1"


rule glnexus:
    input:
        gvcfs=lambda w: expand(
            "results/individual_calls/{sample}.g.vcf.gz",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
    output:
        vcf=temp("results/all_group_samples_joint_calls/{joint_calling_group}.vcf.gz"),
        scratch=temp(
            directory("results/all_group_samples_joint_calls/{joint_calling_group}.DB")
        ),
    params:
        config=config["glnexus"]["config"],
    threads: config["glnexus"]["threads"]
    log:
        "results/logs/glnexus/{joint_calling_group}/stdout.log",
    container:
        "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    shell:
        "glnexus_cli "
        "--config DeepVariantWGS "
        "--dir {output.scratch} "
        "--threads {threads} "
        "{input} "
        "2> {log} "
        "| bcftools view - "
        "| bgzip -c "
        "> {output.vcf} "


rule bcftools_index:
    input:
        "{vcffile}.vcf.gz",
    output:
        "{vcffile}.vcf.gz.csi",
    params:
        extra=config["bcftools_index"]["extra"],
    log:
        "results/logs/bcftools_index/{vcffile}.log",
    threads: config["bcftools_index"]["threads"]
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index "
        "--threads {threads} "
        "--csi "
        "{params.extra} "
        "{input} "
        "1>{log} 2>&1"


rule create_reheader_sample_file:
    input:
        joint_calling_groups=config["joint_calling_groups"],
    output:
        samples=temp("results/joint_calls/{joint_calling_group}_sample_names.tsv"),
    log:
        "results/logs/reheader_sample_file/{joint_calling_group}.log",
    run:
        (
            joint_calling_groups.assign(
                group_sample=lambda x: x.group.str.cat(x.sample_id, sep=":")
            )
            .loc[
                lambda x: x.group == wildcards.joint_calling_group,
                ["sample_id", "group_sample"],
            ]
            .to_csv(output.samples, sep="\t", index=False, header=None)
        )


rule update_sample_names:
    input:
        vcf=rules.glnexus.output.vcf,
        samples=rules.create_reheader_sample_file.output.samples,
    output:
        vcf="results/joint_calls/{joint_calling_group}.vcf.gz",
    log:
        "results/logs/update_sample_names/{joint_calling_group}.log",
    params:
        extra="",
        view_extra="",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "("
        "bcftools reheader "
        "--threads {threads} "
        "--samples {input.samples} "
        "--temp-prefix {input.vcf}.reheader. "
        "{input.vcf} "
        "| bcftools view "
        "--output-type z "
        "--output {output.vcf}; "
        "bcftools index "
        "--csi "
        "{output.vcf} "
        ")1>{log} 2>&1"


rule bcftools_merge:
    input:
        calls=[
            *expand(
                "results/calls/{sample}.vcf.gz",
                sample=(
                    samples.loc[
                        ~samples.sample_id.isin(joint_calling_groups.sample_id)
                    ].sample_id.unique()
                ),
            ),
            *expand(
                "results/individual_calls/{sample}.vcf.gz",
                sample=(
                    samples.loc[
                        samples.sample_id.isin(joint_calling_groups.sample_id)
                    ].sample_id.unique()
                ),
            ),
            *expand(
                "results/joint_calls/{joint_calling_group}.vcf.gz",
                joint_calling_group=joint_calling_group_lists.index,
            ),
        ],
        idxs=[
            *expand(
                "results/calls/{sample}.vcf.gz.csi",
                sample=(
                    samples.loc[
                        ~samples.sample_id.isin(joint_calling_groups.sample_id)
                    ].sample_id.unique()
                ),
            ),
            *expand(
                "results/individual_calls/{sample}.vcf.gz.csi",
                sample=(
                    samples.loc[
                        samples.sample_id.isin(joint_calling_groups.sample_id)
                    ].sample_id.unique()
                ),
            ),
            *expand(
                "results/joint_calls/{joint_calling_group}.vcf.gz.csi",
                joint_calling_group=joint_calling_group_lists.index,
            ),
        ],
    output:
        calls="results/merged_calls/all.unfiltered.vcf.gz",
    log:
        "results/logs/bcftools_merge/bcftools_merge.log",
    params:
        config["bcftools_merge"]["params"],  # optional parameters for bcftools concat (except -o)
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools merge "
        "--output-type z "
        "--output {output.calls} "
        "{params} "
        "{input.calls} "
        "> {log} 2>&1"


rule vcftools_filter:
    input:
        rules.bcftools_merge.output.calls,
    output:
        "results/merged_calls/all.vcf.gz",
    log:
        "results/logs/vcftools_filter.log",
    params:
        filter=config["vcftools_filter"]["filter"],
        extra="--recode --recode-INFO-all",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "(vcftools "
        "--gzvcf {input} "
        "{params.filter} "
        "{params.extra} "
        "--stdout | gzip -c >{output} "
        ")1>{log} 2>&1"


