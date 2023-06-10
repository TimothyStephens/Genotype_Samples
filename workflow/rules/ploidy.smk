

def get_lib_type(wildcards):
	'''
	If any of the units (libraires that are merged into a single bam file) are classified as "rna" then we need to run GATK SplitNCigarReads. 
	Else, we assume they are DNA and we can just use the results directly.
	'''
	use_RNA = False
	rows = samples.loc[(wildcards.sample), ["sample_id", "unit", "lib_type", "fq1", "fq2"]]
	for i, row in rows.iterrows():
		if "rna" in row.lib_type:
			use_RNA = True
	if use_RNA:
		return "rna"
	else:
		return "dna"


rule ploidy_SplitNCigarReads:
	input:
		bam="results/mapping_merged/{ref_name}/{sample}.bam",
		idx="results/mapping_merged/{ref_name}/{sample}.bam.csi",
		ref="resources/{ref_name}/genome.fasta".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
	output:
		bam=temp("results/ploidy/{ref_name}/{sample}.bam"),
		bai=temp("results/ploidy/{ref_name}/{sample}.bai"),
	log:
		"results/logs/ploidy/{ref_name}/SplitNCigarReads/{sample}.log",
	conda:
		"../envs/gatk4.yaml"
	shell:
		"(gatk SplitNCigarReads -R {input.ref} -I {input.bam} -O {output.bam}) 1>{log} 2>&1"


rule ploidy_nQuire_create:
	input:
		bam=rules.ploidy_SplitNCigarReads.input.bam,
		idx=rules.ploidy_SplitNCigarReads.input.idx,
	output:
		"results/ploidy/{ref_name}/{sample}.bin",
	log:
		"results/logs/ploidy/{ref_name}/nQuire_create/{sample}.log",
	params:
		extra=config["ploidy_nQuire"]["create_params"],
		min_quality=config["ploidy_nQuire"]["min_quality"],
		min_coverage=config["ploidy_nQuire"]["min_coverage"],
	container:
		"docker://nanozoo/nquire:0.0--e42aee8"
	shell:
		"(nQuire create -q {params.min_quality} -c {params.min_coverage} -x -b {input.bam} -o $(echo {output} | sed -e 's/.bin//') {params.extra}) 1>{log} 2>&1"


rule ploidy_nQuire_denoise:
	input:
		nQbin=rules.ploidy_nQuire_create.output,
	output:
		"results/ploidy/{ref_name}/{sample}.denoised.bin",
	log:
		"results/logs/ploidy/{ref_name}/nQuire_denoise/{sample}.denoise.log",
	params:
		extra=config["ploidy_nQuire"]["denoise_params"],
	container:
		"docker://nanozoo/nquire:0.0--e42aee8"
	shell:
		"(nQuire denoise -o $(echo {output} | sed -e 's/.bin//') {params.extra} {input.nQbin}) 1>{log} 2>&1"


rule ploidy_nQuire_coverage:
	input:
		bam="results/mapping_merged/{ref_name}/{sample}.bam".format(
			ref_name=list(config["ref_genomes"].keys())[0],
			sample="{sample}",
		),
		nQbin=rules.ploidy_nQuire_denoise.output,
	output:
		"results/ploidy/{ref_name}/{sample}.denoised.bin.coverage.sitesProp.gz",
	log:
		"results/logs/ploidy/{ref_name}/nQuire_denoise/{sample}.nQuire_coverage.log",
	params:
		extra=config["ploidy_nQuire"]["coverage_params"],
	container:
		"docker://nanozoo/nquire:0.0--e42aee8"
	shell:
		"(nQuire view {input.nQbin} -a {input.bam} {params.extra} | awk -F'\\t' '{{print $4/$3\"\\n\"$5/$3}}' | gzip -c > {output}) 1>{log} 2>&1"


rule ploidy_nQuire_site_count:
	input:
		nQbin=rules.ploidy_nQuire_create.output,
		nQbin_denoised=rules.ploidy_nQuire_denoise.output,
	output:
		"results/ploidy/{ref_name}/{sample}.site_counts.tsv",
	log:
		"results/logs/ploidy/{ref_name}/{sample}.count.log",
	container:
		"docker://nanozoo/nquire:0.0--e42aee8"
	shell:
		"("
		"C=$(nQuire view {input.nQbin} | wc -l)"
		"; D=$(nQuire view {input.nQbin_denoised} | wc -l)"
		"; printf \"%s\\t%s\\t%s\\t%.2f\\n\" {wildcards.sample} $C $D $(echo \"scale=2;($D/$C)*100\" | bc)"
		" > {output}"
		")"
		" 1>{log} 2>&1"


rule ploidy_nQuire_merge_site_counts:
	input:
		lambda wildcards: expand("results/ploidy/{ref_name}/{sample}.site_counts.tsv",
			ref_name=list(config["ref_genomes"].keys())[0],
			sample=samples.sample_id.unique()),
	output:
		"results/{project}/ploidy/nQuire_sites_count.txt",
	log:
		"results/logs/{project}/ploidy/nQuire_lrdmodel.log",
	conda:
		"../envs/bash.yaml"
	shell:
		"("
		" echo -e \"sample_id\\tnum_sites_normal\\tnum_sites_denoised\\tprop_denoised\" > {output};"
		" cat {input} >> {output}"
		")"
		" 1>{log} 2>&1"


rule ploidy_nQuire_lrdmodel:
	input:
		nQbins=lambda wildcards: [
			  *expand("results/ploidy/{ref_name}/{sample}.bin", 
				ref_name=list(config["ref_genomes"].keys())[0],
				sample=samples.sample_id.unique()),
			  *expand("results/ploidy/{ref_name}/{sample}.denoised.bin", 
				ref_name=list(config["ref_genomes"].keys())[0],
				sample=samples.sample_id.unique()),
		],
	output:
		"results/{project}/ploidy/nQuire_lrdmodel.txt",
	log:
		"results/logs/{project}/ploidy/nQuire_lrdmodel.log",
	params:
		extra=config["ploidy_nQuire"]["lrdmodel_params"],
	threads: config["ploidy_nQuire"]["threads"]
	container:
		"docker://nanozoo/nquire:0.0--e42aee8"
	shell:
		"(nQuire lrdmodel {params.extra} -t {threads} {input.nQbins} > {output}) 1>{log} 2>&1"


rule format_nQuire_results:
	input:
		rules.ploidy_nQuire_lrdmodel.output,
		rules.ploidy_nQuire_merge_site_counts.output,
	output:
		"results/{project}/final/nQuire.tsv",
	log:
		"results/logs/{project}/ploidy/format_nQuire_results.log",
	params:
		add_values="workflow/scripts/add_value_to_table.py"
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_nQuire_results.sh"


