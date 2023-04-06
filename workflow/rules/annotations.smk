

rule format_annotations:
	input:
		groups=rules.format_results_vcf_clone_detect_groups.output,
		nQuire=rules.format_nQuire_results.output,
	output:
		"results/{project}/final/Annotations.tsv",
	log:
		"results/logs/{project}/relatedness/format_Annotations.log",
	params:
		add_values="workflow/scripts/add_value_to_table.py"
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_Annotations.sh"


rule ploidy_nQuire_overage_file_list:
	input:
		cov=lambda wildcards: expand("results/{project}/ploidy/{sample}.denoised.bin.coverage.sitesProp.gz",
			project=wildcards.project,
			sample=samples.sample_id.unique()),
		annots=rules.format_annotations.output,
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


