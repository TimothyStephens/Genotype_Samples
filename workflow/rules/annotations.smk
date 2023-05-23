

rule format_annotations:
	input:
		samples=config["samples"],
		color_list="workflow/scripts/color_list.txt",
	output:
		samples="results/{project}/final/samples.tsv",
		color_list="results/{project}/final/colors.tsv",
	log:
		"results/logs/{project}/final/format_annotations.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_annotations.sh"


rule combine_genotyping_results:
	input:
		samples=rules.format_annotations.output.samples,
		groups=rules.format_results_vcf_clone_detect_groups.output,
		nQuire=rules.format_nQuire_results.output,
	output:
		samples="results/{project}/final/samples.genotyping.tsv",
	log:
		"results/logs/{project}/relatedness/combine_genotyping_results.log",
	params:
		add_values="workflow/scripts/add_value_to_table.py"
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/combine_genotyping_results.sh"


rule ploidy_nQuire_overage_file_list:
	input:
		cov=lambda wildcards: expand("results/ploidy/{sample}.denoised.bin.coverage.sitesProp.gz",
			sample=samples.sample_id.unique()),
		annots=rules.combine_genotyping_results.output.samples,
	output:
		"results/{project}/ploidy/nQuire_overage_file_list.tsv",
	log:
		"results/logs/{project}/ploidy/nQuire_overage_file_list.log",
	conda:
		"../envs/bash.yaml"
	shell:
		"("
		" awk -F'\\t' 'NR>1' {input.annots} "
		" | sort -k5,5n -k1,1"
		" | awk -F'\\t' '{{"
		"F=$1\".denoised.bin.coverage.sitesProp.gz\";"
		"print F\"\\t\"$1\" (\"$4\")\"}}'"
		" > {output}"
		")"
		" 1>{log} 2>&1"


