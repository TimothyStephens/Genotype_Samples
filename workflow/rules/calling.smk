

rule calling_download_RNAseq_model:
	output:
		"resources/deepvariant_RNAseq_model/model.ckpt"
	log:
		"results/logs/resources/download_RNAseq_model.log"
	shell:
		"("
		"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.data-00000-of-00001 > {output}.data-00000-of-00001;"
		"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.example_info.json > {output}.example_info.json;"
		"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.index > {output}.index;"
		"curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.meta > {output}.meta;"
		"touch {output};"
		") 1>{log} 2>&1"


def get_model_params(wildcards):
	'''
	If any of the units (libraires that are merged into a single bam file) are classified as "dna" then we use the DNA model. 
	Else, we assume they are RNA and we use the RNA model.
	'''
	use_DNA = False
	rows = samples.loc[(wildcards.sample), ["sample_id", "unit", "lib_type", "fq1", "fq2"]]
	for i, row in rows.iterrows():
		if row.lib_type == 'dna':
			use_DNA = True
	if use_DNA:
		return "--model_type WGS"
	else:
		return "--model_type=WES --customized_model="+rules.calling_download_RNAseq_model.output[0]+" --make_examples_extra_args=\"split_skip_reads=true,channels=''\""


rule calling_variants:
	input:
		bam=rules.mapping_merge.output.bam,
		idx=rules.mapping_merge.output.idx,
		ref="resources/{ref_name}/genome.fasta".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
		ref_idx="resources/{ref_name}/genome.fasta.fai".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
		rnaseq_model=rules.calling_download_RNAseq_model.output,
	output:
		vcf="results/{project}/calling/{sample}.vcf.gz",
		gvcf="results/{project}/calling/{sample}.g.vcf.gz",
		report=report(
			"results/{project}/calling/{sample}.visual_report.html",
			caption="../report/vcf.rst",
			subcategory="Deepvariant Reports",
			labels={"Sample ID": "{sample}"},
		),
		tmp_dir=temp(directory("results/{project}/calling/{sample}.tmp")),
	log:
		 "results/{project}/log/calling/{sample}/stdout.log",
	params:
		model_params=lambda wildcards: get_model_params(wildcards),
		extra=config["calling_variants"]["params"],
	threads: config["calling_variants"]["threads"]
	container:
		"docker://google/deepvariant:1.4.0"
	shell:
		"run_deepvariant"
		" {params.model_params}"
		" --ref {input.ref}"
		" --reads {input.bam}"
		" --sample_name {wildcards.sample}"
		" --logging_dir results/{wildcards.project}/log/calling/{wildcards.sample}.deepvariant"
		" --output_vcf {output.vcf}"
		" --output_gvcf {output.gvcf}"
		" --intermediate_results_dir {output.tmp_dir}"
		" --num_shards {threads}"
		" 1>{log} 2>&1"


rule calling_merge_VCFs:
	input:
		calls=lambda wildcards:	expand(
			"results/{project}/calling/{sample}.vcf.gz",
			project=wildcards.project,
			sample=samples.sample_id.unique()
		),
		idxs=lambda wildcards: expand(
			"results/{project}/calling/{sample}.vcf.gz.csi",
			project=wildcards.project,
			sample=samples.sample_id.unique()
		),
	output:
		"results/{project}/merged_calls/calls.unfiltered.vcf.gz",
	log:
		"results/{project}/log/merged_calls/VCFs_merge.log",
	params:
		extra=config["calling_merge_VCFs"]["params"],
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge"
		" --output-type z"
		" --output {output}"
		" --force-samples"
		" {params.extra}"
		" {input.calls}"
		" 1>{log} 2>&1"


rule calling_VCF_index:
	input:
		"results/{project}/{vcffile}.vcf.gz",
	output:
		"results/{project}/{vcffile}.vcf.gz.csi",
	params:
		extra=config["calling_VCFs_index"]["params"],
	log:
		"results/{project}/log/VCF_index/{vcffile}.log",
	threads: config["calling_VCFs_index"]["threads"]
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools index"
		" --threads {threads}"
		" --csi"
		" {params.extra}"
		" {input}"
		" 1>{log} 2>&1"


rule calling_filter_merged_VCF:
	input:
		rules.calling_merge_VCFs.output,
	output:
		"results/{project}/merged_calls/calls.filtered.vcf.gz",
	log:
		"results/{project}/log/merged_calls/VCF_filter.log",
	params:
		filter=config["calling_filter_VCFs"]["filter"],
		extra="--recode --recode-INFO-all",
	conda:
		"../envs/vcftools.yaml"
	shell:
		"("
		"vcftools"
		" --gzvcf {input}"
		" {params.filter}"
		" {params.extra}"
		" --stdout | gzip -c >{output}"
		")"
		" 1>{log} 2>&1"


rule format_filtered_VCF:
	input:
		rules.calling_filter_merged_VCF.output,
	output:
		"results/{project}/final/calls.filtered.vcf.gz"
	log:
		"results/{project}/log/merged_calls/format_filtered_VCF.log",
	shell:
		"(cp {input} {output}) 1>{log} 2>&1"


