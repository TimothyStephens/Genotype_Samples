

rule ploidy_nQuire_create:
	input:
		bam=rules.mapping_merge.output.bam,
		idx=rules.mapping_merge.output.idx,
		programs=rules.install_nQuire.output,
	output:
		"results/{project}/ploidy/{sample}.bin",
	log:
		"results/{project}/log/ploidy/nQuire_create/{sample}.log",
	params:
		extra=config["ploidy_nQuire"]["create_params"],
		min_quality=config["ploidy_nQuire"]["min_quality"],
		min_coverage=config["ploidy_nQuire"]["min_coverage"],
	conda:
		"../envs/nQuire.yaml"
	shell:
		"(nQuire create -q {params.min_quality} -c {params.min_coverage} -x -b {input.bam} -o $(echo {output} | sed -e 's/.bin//') {params.extra}) 1>{log} 2>&1"


rule ploidy_nQuire_denoise:
	input:
		nQbin=rules.ploidy_nQuire_create.output,
		programs=rules.install_nQuire.output,
	output:
		"results/{project}/ploidy/{sample}.denoised.bin",
	log:
		"results/{project}/log/ploidy/nQuire_denoise/{sample}.denoise.log",
	params:
		extra=config["ploidy_nQuire"]["denoise_params"],
	conda:
		"../envs/nQuire.yaml"
	shell:
		"(nQuire denoise -o $(echo {output} | sed -e 's/.bin//') {params.extra} {input.nQbin}) 1>{log} 2>&1"


rule ploidy_nQuire_coverage:
	input:
		bam=rules.mapping_merge.output.bam,
		nQbin=rules.ploidy_nQuire_denoise.output,
		programs=rules.install_nQuire.output,
	output:
		"results/{project}/ploidy/{sample}.denoised.bin.coverage.sitesProp.gz",
	log:
		"results/{project}/log/ploidy/nQuire_denoise/{sample}.nQuire_coverage.log",
	params:
		extra=config["ploidy_nQuire"]["coverage_params"],
	conda:
		"../envs/nQuire.yaml"
	shell:
		"(nQuire view {input.nQbin} -a {input.bam} {params.extra} | awk -F'\\t' '{{print $4/$3\"\\n\"$5/$3}}' | gzip -c > {output}) 1>{log} 2>&1"


rule ploidy_nQuire_site_count:
	input:
		nQbin=rules.ploidy_nQuire_create.output,
		nQbin_denoised=rules.ploidy_nQuire_denoise.output,
		programs=rules.install_nQuire.output,
	output:
		"results/{project}/ploidy/{sample}.site_counts.tsv",
	log:
		"results/{project}/log/ploidy/{sample}.count.log",
	conda:
		"../envs/nQuire.yaml"
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
		lambda wildcards: expand("results/{project}/ploidy/{sample}.site_counts.tsv",
			project=wildcards.project,
			sample=samples.sample_id.unique()),
	output:
		"results/{project}/ploidy/nQuire_sites_count.txt",
	log:
		"results/{project}/log/ploidy/nQuire_lrdmodel.log",
	shell:
		"("
		" echo -e \"sample_id\\tnum_sites_normal\\tnum_sites_denoised\\tprop_denoised\" > {output};"
		" cat {input} >> {output}"
		")"
		" 1>{log} 2>&1"


rule ploidy_nQuire_lrdmodel:
	input:
		nQbins=lambda wildcards: [
			  *expand("results/{project}/ploidy/{sample}.bin", 
				project=wildcards.project,
				sample=samples.sample_id.unique()),
			  *expand("results/{project}/ploidy/{sample}.denoised.bin", 
				project=wildcards.project, 
				sample=samples.sample_id.unique()),
		],
		programs=rules.install_nQuire.output,
	output:
		"results/{project}/ploidy/nQuire_lrdmodel.txt",
	log:
		"results/{project}/log/ploidy/nQuire_lrdmodel.log",
	params:
		extra=config["ploidy_nQuire"]["lrdmodel_params"],
	threads: config["ploidy_nQuire"]["threads"]
	conda:
		"../envs/nQuire.yaml"
	shell:
		"(nQuire lrdmodel {params.extra} -t {threads} {input.nQbins} > {output}) 1>{log} 2>&1"


rule format_nQuire_results:
	input:
		rules.ploidy_nQuire_lrdmodel.output,
		rules.ploidy_nQuire_merge_site_counts.output,
	output:
		"results/{project}/final/nQuire.tsv",
	log:
		"results/{project}/log/ploidy/format_nQuire_results.log",
	params:
		add_values="workflow/scripts/add_value_to_table.py"
	script:
		"../scripts/format_nQuire_results.sh"


