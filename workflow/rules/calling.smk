

rule deepvariant:
    input:
        bam=rules.samtools_merge.output,
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
        top="results/logs/deepvariant/{sample}/stdout.log",
        make_examples="results/logs/deepvariant/{sample}/make_examples.log",
	call_variants="results/logs/deepvariant/{sample}/call_variants.log",
        postprocess_variants="results/logs/deepvariant/{sample}/postprocess_variants.log",
    conda:
        "../envs/deepvariant.yaml"
    shell:
        "( rm -fr {params.tmp_dir}; "
        "mkdir -p {params.tmp_dir}; "
        "dv_make_examples.py "
        "--cores {threads} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "--sample {wildcards.sample} "
        "--examples {params.tmp_dir} "
        "--logdir results/logs/deepvariant/{wildcards.sample} "
        "{params.extra} "
        "1>{log.make_examples} 2>&1; "
        "dv_call_variants.py "
        "--cores {threads} "
        "--outfile {params.tmp_dir}/{wildcards.sample}.tmp "
        "--sample {wildcards.sample} "
        "--examples {params.tmp_dir} "
        "--model {params.model} "
        " 1>{log.call_variants} 2>&1; "
        "dv_postprocess_variants.py "
        "--ref {input.ref} "
        "--infile {params.tmp_dir}/{wildcards.sample}.tmp "
        "--outfile {output.vcf} "
        "1>{log.postprocess_variants} 2>&1"
        ") 1>{log.top} 2>&1"


rule deepvariant_gvcf:
    input:
        bam=rules.samtools_merge.output,
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
        top="results/logs/deepvariant_gvcf/{sample}/stdout.log",
        make_examples="results/logs/deepvariant_gvcf/{sample}/make_examples.log",
        call_variants="results/logs/deepvariant_gvcf/{sample}/call_variants.log",
        postprocess_variants="results/logs/deepvarianti_gvcf/{sample}/postprocess_variants.log",
    conda:
        "../envs/deepvariant.yaml"
    shell:
        "( rm -fr {params.tmp_dir}; "
        "mkdir -p {params.tmp_dir}; "
        "dv_make_examples.py "
        "--cores {threads} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "--sample {wildcards.sample} "
        "--examples {params.tmp_dir} "
        "--logdir results/logs/deepvariant_gvcf/{wildcards.sample} "
        "--gvcf {params.tmp_dir} "
        "{params.extra} "
        "1>{log.make_examples} 2>&1; "
        "dv_call_variants.py "
        "--cores {threads} "
        "--outfile {params.tmp_dir}/{wildcards.sample}.tmp "
        "--sample {wildcards.sample} "
        "--examples {params.tmp_dir} "
        "--model {params.model} "
        " 1>{log.call_variants} 2>&1; "
        "dv_postprocess_variants.py "
        "--ref {input.ref} "
        "--gvcf_infile {params.tmp_dir}/{wildcards.sample}.gvcf.tfrecord@{threads}.gz "
        "--gvcf_outfile {output.gvcf} "
        "--infile {params.tmp_dir}/{wildcards.sample}.tmp "
        "--outfile {output.vcf} "
        "1>{log.postprocess_variants} 2>&1"
        ") 1>{log.top} 2>&1"


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
        extra=config["bcftools_index"]["extra"] + " --threads {}".format(
            config["bcftools_index"]["threads"]
        ),
    log:
        "results/logs/bcftools_index/{vcffile}.log",
    threads: config["bcftools_index"]["threads"]
    wrapper:
        "v1.23.4/bio/bcftools/index"


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
        view_extra="-O z",
    wrapper:
        "v1.23.4/bio/bcftools/reheader"


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
        config["bcftools_merge"]["params"] + " -Oz",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "v1.23.4/bio/bcftools/merge"


rule bcftools_filter:
    input:
        rules.bcftools_merge.output.calls,
    output:
        "results/merged_calls/all.vcf.gz",
    log:
        "results/logs/bcftools_filter_all.log",
    params:
        filter=config["bcftools_filter"]["filter"],
        extra="",
    wrapper:
        "v1.23.4/bio/bcftools/filter"


