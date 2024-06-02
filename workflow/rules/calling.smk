

rule calling_download_RNAseq_model:
	output:
		"resources/deepvariant_RNAseq_model/model.ckpt"
	log:
		"results/logs/resources/download_RNAseq_model.log"
	conda:
		"../envs/bash.yaml"
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
		if 'dna' in row.lib_type:
			use_DNA = True
	if use_DNA:
		return "--model_type WGS"
	else:
		return "--model_type=WES --customized_model="+rules.calling_download_RNAseq_model.output[0]+" --make_examples_extra_args=\"split_skip_reads=true,channels=''\""


rule calling_variants:
	input:
		bam=rules.mapping_merge.output.bam,
		idx=rules.mapping_merge.output.idx,
		ref=rules.ref_parse_genome.output,
		ref_idx=rules.ref_faidx.output,
		rnaseq_model=rules.calling_download_RNAseq_model.output,
	output:
		vcf="results/calling/{ref_name}/{sample}.vcf.gz",
		gvcf="results/calling/{ref_name}/{sample}.g.vcf.gz",
		report=report(
			"results/calling/{ref_name}/{sample}.visual_report.html",
			caption="../report/vcf.rst",
			subcategory="Deepvariant Reports",
			labels={"Sample ID": "{sample}"},
		),
		tmp_dir=temp(directory("results/calling/{ref_name}/{sample}.tmp")),
	log:
		 "results/logs/calling/{ref_name}/{sample}/stdout.log",
	params:
		model_params=lambda wildcards: get_model_params(wildcards),
		extra=config["calling_variants"]["params"],
	threads: config["calling_variants"]["threads"]
	resources:
		mem_gb=config["calling_variants"]["memory"]
	container:
		"docker://google/deepvariant:1.4.0"
	shell:
		"run_deepvariant"
		" {params.model_params}"
		" --ref {input.ref}"
		" --reads {input.bam}"
		" --sample_name {wildcards.sample}"
		" --logging_dir results/logs/calling/{wildcards.ref_name}/{wildcards.sample}.deepvariant"
		" --output_vcf {output.vcf}"
		" --output_gvcf {output.gvcf}"
		" --intermediate_results_dir {output.tmp_dir}"
		" --num_shards {threads}"
		" 1>{log} 2>&1"


rule calling_merge_VCFs:
	input:
		calls=lambda wildcards:	expand(
			"results/calling/{ref_name}/{sample}.g.vcf.gz",
			ref_name=list(config["ref_genomes"].keys())[0],
			sample=samples.sample_id.unique()
		),
		idxs=lambda wildcards: expand(
			"results/calling/{ref_name}/{sample}.g.vcf.gz.csi",
			ref_name=list(config["ref_genomes"].keys())[0],
			sample=samples.sample_id.unique()
		),
	output:
		"results/{project}/calling_merged/calls.unfiltered.vcf.gz",
	log:
		"results/logs/{project}/calling_merged/VCFs_merge.log",
	params:
		extra=config["calling_merge_VCFs"]["params"],
		mem_gb=config["calling_merge_VCFs"]["memory"],
		tmp=temp(directory("results/{project}/calling_merged/GLnexus.DB")),
		bcf=temp("results/{project}/calling_merged/calls.unfiltered.bcf"),
	threads: config["calling_merge_VCFs"]["threads"]
	container:
		"docker://quay.io/mlin/glnexus:v1.3.1"
	shell:
		"("
		"glnexus_cli"
		" {params.extra}"
		" --config DeepVariant"
		" --dir {params.tmp}"
		" --mem-gbytes {params.mem_gb}"
		" --threads {threads}"
		" {input.calls}"
		" > {params.bcf};"
		"bcftools view {params.bcf} | gzip -c > {output}"
		") 1>{log} 2>&1"


rule calling_merged_VCF_stats:
	input:
		rules.calling_merge_VCFs.output,
	output:
		lqual  = "results/{project}/calling_merged/calls.unfiltered.lqual",
		ldepth = "results/{project}/calling_merged/calls.unfiltered.ldepth.mean",
		lmiss  = "results/{project}/calling_merged/calls.unfiltered.lmiss",
		frq    = "results/{project}/calling_merged/calls.unfiltered.frq",
		idepth = "results/{project}/calling_merged/calls.unfiltered.idepth",
		imiss  = "results/{project}/calling_merged/calls.unfiltered.imiss",
		het    = "results/{project}/calling_merged/calls.unfiltered.het",
	log:
		"results/logs/{project}/calling_merged/merged_VCF_stats.log",
	params:
		out="results/{project}/calling_merged/calls.unfiltered",
	conda:
		"../envs/vcftools.yaml"
	shell:
		"("
		# Calculate allele frequency
		"vcftools --gzvcf {input} --freq2 --out {params.out} --max-alleles 2; "
		
		# Calculate mean depth per individual
		"vcftools --gzvcf {input} --depth --out {params.out}; "
		
		# Calculate mean depth per site
		"vcftools --gzvcf {input} --site-mean-depth --out {params.out}; "
		
		# Calculate site quality
		"vcftools --gzvcf {input} --site-quality --out {params.out}; "
		
		# Calculate proportion of missing data per individual
		"vcftools --gzvcf {input} --missing-indv --out {params.out}; "
		
		# Calculate proportion of missing data per site
		"vcftools --gzvcf {input} --missing-site --out {params.out}; "
		
		# Calculate heterozygosity and inbreeding coefficient per individual
		"vcftools --gzvcf {input} --het --out {params.out} "
		
		") 1>{log} 2>&1"


rule calling_VCF_index:
	input:
		"results/{project}/{vcffile}.vcf.gz",
	output:
		"results/{project}/{vcffile}.vcf.gz.csi",
	params:
		extra=config["calling_VCFs_index"]["params"],
	log:
		"results/logs/{project}/VCF_index/{vcffile}.log",
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
		"results/{project}/calling_merged/calls.filtered.vcf.gz",
	log:
		"results/logs/{project}/calling_merged/VCF_filter.log",
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


rule calling_merged_filtered_VCF_stats:
	input:
		rules.calling_filter_merged_VCF.output,
	output:
		lqual  = "results/{project}/calling_merged/calls.filtered.lqual",
		ldepth = "results/{project}/calling_merged/calls.filtered.ldepth.mean",
		lmiss  = "results/{project}/calling_merged/calls.filtered.lmiss",
		frq    = "results/{project}/calling_merged/calls.filtered.frq",
		idepth = "results/{project}/calling_merged/calls.filtered.idepth",
		imiss  = "results/{project}/calling_merged/calls.filtered.imiss",
		het    = "results/{project}/calling_merged/calls.filtered.het",
	log:
		"results/logs/{project}/calling_merged/filtered_VCF_stats.log",
	params:
		out="results/{project}/calling_merged/calls.filtered",
	conda:
		"../envs/vcftools.yaml"
	shell:
		"("
		# Calculate allele frequency
		"vcftools --gzvcf {input} --freq2 --out {params.out} --max-alleles 2; "
		
		# Calculate mean depth per individual
		"vcftools --gzvcf {input} --depth --out {params.out}; "
		
		# Calculate mean depth per site
		"vcftools --gzvcf {input} --site-mean-depth --out {params.out}; "
		
		# Calculate site quality
		"vcftools --gzvcf {input} --site-quality --out {params.out}; "
		
		# Calculate proportion of missing data per individual
		"vcftools --gzvcf {input} --missing-indv --out {params.out}; "
		
		# Calculate proportion of missing data per site
		"vcftools --gzvcf {input} --missing-site --out {params.out}; "
		
		# Calculate heterozygosity and inbreeding coefficient per individual
		"vcftools --gzvcf {input} --het --out {params.out} "
		
		") 1>{log} 2>&1"


rule format_filtered_VCF:
	input:
		rules.calling_filter_merged_VCF.output,
	output:
		"results/{project}/final/calls.filtered.vcf.gz"
	log:
		"results/logs/{project}/calling_merged/format_filtered_VCF.log",
	conda:
		"../envs/bash.yaml"
	shell:
		"(cp {input} {output}) 1>{log} 2>&1"


