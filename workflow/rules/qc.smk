

rule qc_fastqc:
	input:
		"results/trimmed/{fq}.fastq.gz",
	output:
		html="results/qc/trimmed/{fq}_fastqc.html",
		zip="results/qc/trimmed/{fq}_fastqc.zip",
		tmpdir=temp(directory("results/qc/trimmed/{fq}.fastqc_tmpdir")),
	log:
		"results/logs/qc/trimmed/{fq}.log",
	params:
		extra=config["qc_FastQC"]["params"],
	threads: config["qc_FastQC"]["threads"]
	conda:
		"../envs/fastqc.yaml"
	shell:
		"("
		"prefix=$(basename {input} | sed -e 's/.fastq.*//' -e 's/.fq.*//'); "
		"rm -fr {output.tmpdir}; mkdir -p {output.tmpdir}; "
		"fastqc"
		" {params.extra}"
		" -t {threads}"
		" --outdir {output.tmpdir}"
		" {input};"
		" mv {output.tmpdir}/${{prefix}}_fastqc.zip  {output.zip};"
		" mv {output.tmpdir}/${{prefix}}_fastqc.html {output.html};"
		")"
		" 1>{log} 2>&1"


def get_raw_fastq_paths(wildcards):
	for i, row in samples.iterrows():
		if wildcards.fq == "{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type):
			return [row.fq1]
		else:
			return [row.fq2]


rule qc_fastqc_rawReads:
	input:
		unpack(get_raw_fastq_paths),
	output:
		html="results/qc/raw/{fq}_fastqc.html",
		zip="results/qc/raw/{fq}_fastqc.zip",
		tmpdir=temp(directory("results/qc/raw/{fq}.fastqc_tmpdir")),
	log:
		"results/logs/qc/raw/{fq}.log",
	params:
		extra=config["qc_FastQC"]["params"],
	threads: config["qc_FastQC"]["threads"]
	conda:
		"../envs/fastqc.yaml"
	shell:
		"("
		" in_fq=$(basename {input});"
		" in_prefix=$(echo ${{in_fq}} | sed -e 's/.fq.*//' -e 's/.fastq.*//');"
		" out_fq=$(basename {output.html} | sed -e 's/.html/.fastq.gz/');"
		" out_prefix=$(echo ${{out_fq}} | sed -e 's/.fq.*//' -e 's/.fastq.*//');"
		" echo in_fq=${{in_fq}} in_prefix=${{in_prefix}} out_fq=${{out_fq}} out_prefix=${{out_prefix}};"
		""
		" rm -fr {output.tmpdir};"
		" mkdir -p {output.tmpdir};"
		""
		"fastqc"
		" {params.extra}"
		" -t {threads}"
		" --extract"
		" --outdir {output.tmpdir}"
		" {input};"
		""
		" ("
		" cd {output.tmpdir}/;"
		" sed -i -e \"s/${{in_fq}}/${{out_fq}}/\" ${{in_prefix}}_fastqc/*.txt ${{in_prefix}}_fastqc/*.fo ${{in_prefix}}_fastqc/*.html ${{in_prefix}}_fastqc.html;"
		" mv ${{in_prefix}}_fastqc ${{out_prefix}}_fastqc;"
		" mv ${{in_prefix}}_fastqc.html ${{out_prefix}}_fastqc.html;"
		" zip -r ${{out_prefix}}_fastqc.zip ${{out_prefix}}_fastqc;"
		" );"
		""
		" mv {output.tmpdir}/${{out_prefix}}_fastqc.zip  {output.zip};"
		" mv {output.tmpdir}/${{out_prefix}}_fastqc.html {output.html};"
		")"
		" 1>{log} 2>&1"


rule qc_samtools_stats:
	input:
		bam="results/"+PROJECT+"/mapping_merged/{sample}.bam",
	output:
		"results/"+PROJECT+"/qc/mapping_merged/{sample}_samtools_stats.txt",
	log:
		"results/"+PROJECT+"/log/qc/mapping_merged/{sample}.log",
	params:
		extra=config["qc_samtools_stats"]["params"],
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools stats"
		" {params.extra}"
		" {input.bam}"
		" 1>{output}"
		" 2>{log}"


rule qc_qualimap:
	input:
		bam="results/"+PROJECT+"/mapping_merged/{sample}.bam",
	output:
		directory("results/"+PROJECT+"/qc/mapping_merged/{sample}_qualimap"),
	log:
		"results/"+PROJECT+"/log/qc/mapping_merged/{sample}.qualimap.log",
	params:
		extra=config["qc_qualimap"]["params"],
	conda:
		"../envs/qualimap.yaml"
	shell:
		"qualimap bamqc"
		" {params.extra}"
		" -bam {input.bam}"
		" -outdir {output}"
		" -outformat HTML"
		" 1>{log} 2>&1"


rule qc_samtools_stats_unmerged:
	input:
		bam="results/"+PROJECT+"/mapping/{fq}.sorted.bam",
	output:
		"results/"+PROJECT+"/qc/mapping/{fq}_samtools_stats.txt",
	log:
		"results/"+PROJECT+"/log/qc/mapping/{fq}.log",
	params:
		extra=config["qc_samtools_stats"]["params"],
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools stats"
		" {params.extra}"
		" {input.bam}"
		" 1>{output}"
		" 2>{log}"


rule qc_bcftools_stats:
	input:
		"results/"+PROJECT+"/merged_calls/{name}.vcf.gz",
	output:
		"results/"+PROJECT+"/qc/merged_calls/{name}_bcftools_stats.txt",
	log:
		"results/"+PROJECT+"/log/qc/merged_calls/{name}_bcftools_stats.log",
	params:
		extra=config["qc_bcftools_stats"]["params"],
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools stats"
		" {params.extra}"
		" {input}"
		" > {output}"
		" 2>{log}"


def expand_fastq_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2):
			out.append("{lib_type}/pe/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("{lib_type}/pe/{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		else:
			out.append("{lib_type}/se/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


def expand_sample_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2):
			out.append("{lib_type}/pe/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		else:
			out.append("{lib_type}/se/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out      


def expand_raw_fastqc_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2):
			out.append("{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if pd.isnull(row.fq2):
			out.append("{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out



rule qc_multiqc_reads:
	input:
		expand("results/qc/raw/{fq}_fastqc.zip", fq=expand_raw_fastqc_paths()),
		expand("results/qc/trimmed/{fq}_fastqc.zip", fq=expand_fastq_paths()),
		expand("results/qc/trimmed/{fq}_fastp.json", fq=expand_sample_paths()),
		expand("results/"+PROJECT+"/qc/mapping/{fq}_samtools_stats.txt", fq=expand_sample_paths()),
		#expand("results/"+PROJECT+"/qc/mapping_merged/{sample}_samtools_stats.txt", sample=samples.sample_id.unique()),
		expand("results/"+PROJECT+"/qc/mapping_merged/{sample}_qualimap", sample=samples.sample_id.unique()),
	output:
		report(
			"results/"+PROJECT+"/qc/multiqc_reads.html",
			caption="../report/multiqc_reads.rst",
			category=PROJECT,
			subcategory="MultiQC",
			labels={"QC": "Reads"},
		),
	log:
		"results/"+PROJECT+"/log/qc/multiqc_reads.log",
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


rule qc_multiqc_variant_calls:
	input:
		expand("results/"+PROJECT+"/qc/merged_calls/{name}_bcftools_stats.txt", name=["calls.unfiltered", "calls.filtered"]),
	output:
		report(
			"results/"+PROJECT+"/qc/multiqc_calls.html",
			caption="../report/multiqc_calls.rst",
			category=PROJECT,
			subcategory="MultiQC",
			labels={"QC": "Called variants"},
		),
	log:
		"results/"+PROJECT+"/log/qc/multiqc_calls.log",
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


