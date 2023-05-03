

rule multiqc_reads:
	input:
		expand("results/qc/raw/{fq}_fastqc.zip", fq=expand_raw_fastqc_paths()),
		expand("results/qc/trimmed/{fq}_fastqc.zip", fq=expand_fastq_paths()),
		expand("results/qc/trimmed/{fq}_fastp.json", fq=expand_sample_paths()),
		lambda wildcards: expand("results/qc/mapping/{ref_name}/{fq}_samtools_stats.txt", 
			ref_name=list(config["ref_genomes"].keys())[0],
			fq=expand_sample_paths()
		),
		lambda wildcards: expand("results/qc/mapping_merged/{ref_name}/{sample}_qualimap", 
			ref_name=list(config["ref_genomes"].keys())[0], 
			sample=samples.sample_id.unique()
		),
	output:
		report(
			"results/{project}/qc/multiqc/reads.html",
			caption="../report/multiqc_reads.rst",
			subcategory="MultiQC",
			labels={"QC": "Reads"},
		),
	log:
		"results/logs/{project}/qc/multiqc/reads.log",
	params:
		extra=config["qc_multiqc_reads"]["params"],
	conda:
		"../envs/multiqc.yaml"
	shell:
		"output_dir=$(dirname {output}); "
		"output_name=$(basename {output}); "
		"multiqc"
		" {params.extra}"
		" --config workflow/report/multiqc_reads_config.yaml"
		" --force"
		" -o $output_dir"
		" -n $output_name"
		" {input}"
		" 1>{log} 2>&1"


rule multiqc_variant_calls:
	input:
		lambda wildcards: expand("results/{project}/qc/calling_merged/{name}_bcftools_stats.txt", 
			project=wildcards.project, 
			name=["calls.unfiltered", "calls.filtered"]),
	output:
		report(
			"results/{project}/qc/multiqc/calls.html",
			caption="../report/multiqc_calls.rst",
			subcategory="MultiQC",
			labels={"QC": "Called variants"},
		),
	log:
		"results/logs/{project}/qc/multiqc/calls.log",
	params:
		extra=config["qc_multiqc_calls"]["params"],
	conda:
		"../envs/multiqc.yaml"
	shell:
		"output_dir=$(dirname {output}); "
		"output_name=$(basename {output}); "
		"multiqc"
		" {params.extra}"
		" --config workflow/report/multiqc_calls_config.yaml"
		" --force"
		" -o $output_dir"
		" -n $output_name"
		" {input}"
		" 1>{log} 2>&1"


rule multiqc_crossMapping:
	input:
		expand("results/qc/raw/{fq}_fastqc.zip", fq=expand_raw_fastqc_paths()),
		expand("results/qc/trimmed/{fq}_fastqc.zip", fq=expand_fastq_paths()),
		expand("results/qc/trimmed/{fq}_fastp.json", fq=expand_sample_paths()),
		lambda wildcards: expand("results/cross-mapping/{s}.samtools_stats.txt",
			s=expand_crossMapping_results_paths()),
	output:
		report(
			"results/{project}/qc/multiqc_crossMapping.html",
			caption="../report/multiqc_crossMapping.rst",
			subcategory="MultiQC",
			labels={"QC": "Reads"},
		),
	log:
		"results/logs/{project}/qc/multiqc/crossMapping_reads.log",
	params:
		extra=config["qc_multiqc_crossMapping"]["params"],
	conda:
		"../envs/multiqc.yaml"
	shell:
		"output_dir=$(dirname {output}); "
		"output_name=$(basename {output}); "
		"multiqc"
		" {params.extra}"
		" --config workflow/report/multiqc_crossMapping_config.yaml"
		" --force"
		" -o $output_dir"
		" -n $output_name"
		" {input}"
		" 1>{log} 2>&1"


