

rule format_annotations:
	input:
		groups=rules.format_results_vcf_clone_detect_groups.output,
		nQuire=rules.format_nQuire_results.output,
	output:
		"results/"+PROJECT+"/final/Annotations.tsv",
	log:
		"results/"+PROJECT+"/log/relatedness/format_Annotations.log",
	params:
		add_values="workflow/scripts/add_value_to_table.py"
	script:
		"../scripts/format_Annotations.sh"


