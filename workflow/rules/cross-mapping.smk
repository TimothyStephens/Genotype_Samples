

rule crossMapping_DNA_pe:
	input:
		reads=rules.trimming_DNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta",
		idx_build=multiext("resources/{ref_name}/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		"results/"+PROJECT+"/cross-mapping/{ref_name}/dna/pe/{sample}-{unit}.samtools_stats.txt"
	log:
		"results/"+PROJECT+"/log/cross-mapping/{ref_name}/dna/pe/{sample}-{unit}.log",
	params:
		mapping_extra=config["crossMapping_DNA_pe"]["mapping_params"],
		sort_extra=config["crossMapping_DNA_pe"]["sort_params"],
		stats_extra=config["crossMapping_DNA_pe"]["stats_params"],
		tmpdir=temp(directory("results/"+PROJECT+"/cross-mapping/{ref_name}/dna/pe/{sample}-{unit}.samtools_tmp")),
	threads: config["crossMapping_DNA_pe"]["threads"]
	conda:
		"../envs/bwa-mem2.yaml"
	shell:
		"("
		"bwa-mem2 mem"
		" -t {threads}"
		" {params.mapping_extra}"
		" {input.idx}"
		" {input.reads}"
		" | samtools sort"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
		" | samtools stats"
		" {params.stats_extra}"
		" -"
		" 1>{output}"
		" && rm -fr {params.tmpdir}"
		")"
		" 1>{log} 2>&1"


rule crossMapping_DNA_se:
	input:
		reads=rules.trimming_DNA_se.output.trimmed,
		idx="resources/{ref_name}/genome.fasta",
		idx_build=multiext("resources/{ref_name}/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		"results/"+PROJECT+"/cross-mapping/{ref_name}/dna/se/{sample}-{unit}.samtools_stats.txt"
	log:
		"results/"+PROJECT+"/log/cross-mapping/{ref_name}/dna/se/{sample}-{unit}.log",
	params:
		mapping_extra=config["crossMapping_DNA_se"]["mapping_params"],
		sort_extra=config["crossMapping_DNA_se"]["sort_params"],
		stats_extra=config["crossMapping_DNA_se"]["stats_params"],
		tmpdir=temp(directory("results/"+PROJECT+"/cross-mapping/{ref_name}/dna/se/{sample}-{unit}.samtools_tmp")),
	threads: config["crossMapping_DNA_se"]["threads"]
	conda:
		"../envs/bwa-mem2.yaml"
	shell:
		"("
		"bwa-mem2 mem"
		" -t {threads}"
		" {params.mapping_extra}"
		" {input.idx}"
		" {input.reads}"
		" | samtools sort"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
		" | samtools stats"
		" {params.stats_extra}"
		" -"
		" 1>{output}"
		" && rm -fr {params.tmpdir}"
		")"
		" 1>{log} 2>&1"


rule crossMapping_RNA_pe:
	input:
		reads=rules.trimming_RNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta.STAR",
	output:
		stats="results/"+PROJECT+"/cross-mapping/{ref_name}/rna/pe/{sample}-{unit}.samtools_stats.txt",
		tmpdir=temp(directory("results/"+PROJECT+"/cross-mapping/{ref_name}/rna/pe/{sample}-{unit}.STAR")),
	log:
		"results/"+PROJECT+"/log/cross-mapping/{ref_name}/rna/pe/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["crossMapping_RNA_pe"]["sjdbOverhang"],
		mapping_extra=config["crossMapping_RNA_pe"]["mapping_params"],
		stats_extra=config["crossMapping_RNA_pe"]["stats_params"],
	threads: config["crossMapping_RNA_pe"]["threads"]
	conda:
		"../envs/star.yaml"
	shell:
		"("
		"rm -fr {output.tmpdir}/tmp; "
		"ulimit -n 100000; "
		"STAR"
		" --runThreadN {threads}"
		" --genomeDir {input.idx}"
		" --readFilesIn {input.reads}"
		" --readFilesCommand \"gunzip -c\""
		" --outSAMtype BAM SortedByCoordinate"
		" --outSAMunmapped Within"
		" --sjdbOverhang {params.sjdbOverhang}"
		" --outFileNamePrefix {output.tmpdir}/STAR."
		" {params.mapping_extra}"
		" --outTmpDir {output.tmpdir}/tmp"
		" --outStd BAM_SortedByCoordinate"
		" | samtools stats"
		" {params.stats_extra}"
		" -"
		" 1>{output.stats}"
		")"
		" 1>{log} 2>&1"


rule crossMapping_RNA_se:
	input:
		reads=rules.trimming_RNA_se.output.trimmed,
		idx="resources/{ref_name}/genome.fasta.STAR",
	output:
		stats="results/"+PROJECT+"/cross-mapping/{ref_name}/rna/se/{sample}-{unit}.samtools_stats.txt",
		tmpdir=temp(directory("results/"+PROJECT+"/cross-mapping/{ref_name}/rna/se/{sample}-{unit}.STAR")),
	log:
		"results/"+PROJECT+"/log/cross-mapping/{ref_name}/rna/se/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["crossMapping_RNA_se"]["sjdbOverhang"],
		mapping_extra=config["crossMapping_RNA_se"]["mapping_params"],
		stats_extra=config["crossMapping_RNA_se"]["stats_params"],
	threads: config["crossMapping_RNA_se"]["threads"]
	conda:
		"../envs/star.yaml"
	shell:
		"("
		"rm -fr {output.tmpdir}/tmp; "
		"ulimit -n 100000; "
		"STAR"
		" --runThreadN {threads}"
		" --genomeDir {input.idx}"
		" --readFilesIn {input.reads}"
		" --readFilesCommand \"gunzip -c\""
		" --outSAMtype BAM SortedByCoordinate"
		" --outSAMunmapped Within"
		" --sjdbOverhang {params.sjdbOverhang}"
		" --outFileNamePrefix {output.tmpdir}/STAR."
		" {params.mapping_extra}"
		" --outTmpDir {output.tmpdir}/tmp"
		" --outStd BAM_SortedByCoordinate"
		" | samtools stats"
		" {params.stats_extra}"
		" -"
		" 1>{output.stats}"
		")"
		" 1>{log} 2>&1"


def expand_crossMapping_results_paths():
	out = []
	for ref_name in list(config["ref_genomes"].keys()):
		for i, row in samples.iterrows():
			if pd.notnull(row.fq2):
				out.append("{ref_name}/{lib_type}/pe/{sample}-{unit}".format(
					sample=row.sample_id, 
					unit=row.unit, 
					lib_type=row.lib_type,
					project=PROJECT,
					ref_name=ref_name,
				))
			else:
				out.append("{ref_name}/{lib_type}/se/{sample}-{unit}".format(
						sample=row.sample_id, 
						unit=row.unit, 
						lib_type=row.lib_type,
						project=PROJECT,
						ref_name=ref_name,
					))
	return out


rule format_crossMapping_results:
	input:
		expand("results/"+PROJECT+"/cross-mapping/{s}.samtools_stats.txt", s=expand_crossMapping_results_paths()),
	output:
		"results/"+PROJECT+"/final/mapping_rates.tsv"
	log:
		"results/"+PROJECT+"/log/combine_mapping_results.log",
	script:
		"../scripts/format_crossMapping_results.sh"


rule qc_multiqc_crossMapping:
	input:
		expand("results/qc/raw/{fq}_fastqc.zip", fq=expand_raw_fastqc_paths()),
		expand("results/qc/trimmed/{fq}_fastqc.zip", fq=expand_fastq_paths()),
		expand("results/qc/trimmed/{fq}_fastp.json", fq=expand_sample_paths()),
		expand("results/"+PROJECT+"/cross-mapping/{s}.samtools_stats.txt", s=expand_crossMapping_results_paths()),
	output:
		report(
			"results/"+PROJECT+"/qc/multiqc_crossMapping.html",
			caption="../report/qc/multiqc_crossMapping.rst",
			category=PROJECT,
			subcategory="MultiQC",
			labels={"QC": "Reads"},
		),
	log:
		"results/"+PROJECT+"/log/qc/multiqc_crossMapping.log",
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


rule plotting_cross_mapping_results:
	input:
		matrix=rules.format_crossMapping_results.output,
	output:
		report(
			"results/"+PROJECT+"/final/cross_mapping_results.html",
			#caption="../report/multiqc_calls.rst",
			category=PROJECT,
			subcategory="Result Plots",
			labels={"Results": "Cross mapping results"},
		),
	log:
		"results/"+PROJECT+"/log/relatedness/plot_cross_mapping_results.log",
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/plot_cross_mapping_results.Rmd"

