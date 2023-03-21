

rule qc_fastqc:
	input:
		"results/trimmed/{fq}.fastq.gz",
	output:
		html="results/qc/fastqc/{fq}_fastqc.html",
		zip="results/qc/fastqc/{fq}_fastqc.zip",
		tmpdir=temp(directory("results/qc/fastqc/{fq}.fastqc_tmpdir")),
	log:
		"results/logs/fastqc/{fq}.log",
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
		html="results/qc/fastqc_raw/{fq}_fastqc.html",
		zip="results/qc/fastqc_raw/{fq}_fastqc.zip",
		tmpdir=temp(directory("results/qc/fastqc_raw/{fq}.fastqc_tmpdir")),
	log:
		"results/logs/fastqc_raw/{fq}.log",
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
		bam="results/mapped/{sample}.bam",
	output:
		"results/qc/samtools_stats/{sample}_bam_stats_merged.txt",
	log:
		"results/logs/samtools_stats/{sample}.log",
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
		bam="results/mapped/{sample}.bam",
	output:
		directory("results/qc/qualimap/{sample}"),
	log:
		"results/logs/qualimap/{sample}.log",
	params:
		extra=config["qc_qualimap"]["params"],
	conda:
		"../envs/qualimap.yaml"
	shell:
		"qualimap bamqc"
		" {params.extra}"
		" -bam {input.bam}"
		" -outdir {output}"
		" 1>{log} 2>&1"


rule qc_samtools_stats_unmerged:
	input:
		bam="results/mapped/{fq}.sorted.bam",
	output:
		"results/qc/samtools_stats/{fq}_bam_stats_unmerged.txt",
	log:
		"results/logs/samtools_stats/{fq}.log",
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
		"results/merged_calls/{name}.vcf.gz",
	output:
		"results/qc/bcftools_stats/{name}_stats.txt",
	log:
		"results/logs/bcftools_stats/{name}.log",
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
		if row.lib_type == "dna" and pd.notnull(row.fq2):
			out.append("dna/pe/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("dna/pe/{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "dna" and pd.isnull(row.fq2):
			out.append("dna/se/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "rna" and pd.notnull(row.fq2):
			out.append("rna/pe/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("rna/pe/{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "rna" and pd.isnull(row.fq2):
			out.append("rna/se/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


def expand_sample_paths():
	out = []
	for i, row in samples.iterrows():
		if row.lib_type == "dna" and pd.notnull(row.fq2):
			out.append("dna/pe/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "dna" and pd.isnull(row.fq2):
			out.append("dna/se/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "rna" and pd.notnull(row.fq2):
			out.append("rna/pe/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "rna" and pd.isnull(row.fq2):
			out.append("rna/se/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
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



rule qc_multiqc:
	input:
		expand("results/qc/fastqc/{fq}_fastqc.zip", fq=expand_fastq_paths()),
		expand("results/qc/fastqc_raw/{fq}_fastqc.zip", fq=expand_raw_fastqc_paths()),
		expand("results/qc/trimmed/{fq}_fastp.json", fq=expand_sample_paths()),
		expand("results/qc/samtools_stats/{sample}_bam_stats_merged.txt", sample=samples.sample_id.unique()),
		expand("results/qc/samtools_stats/{fq}_bam_stats_unmerged.txt", fq=expand_sample_paths()),
		expand("results/qc/qualimap/{sample}", sample=samples.sample_id.unique()),
		expand("results/qc/bcftools_stats/{name}_stats.txt", name=["all.unfiltered", "all"]),
	output:
		report(
			"results/qc/multiqc.html",
			caption="../report/multiqc.rst",
			category="Quality control",
		),
	log:
		"results/logs/multiqc.log",
	params:
		extra=config["qc_multiqc"]["params"],
	conda:
		"../envs/multiqc.yaml"
	shell:
		"output_dir=$(dirname {output}); "
		"output_name=$(basename {output}); "
		"multiqc"
		" {params.extra}"
		" --config workflow/report/multiqc_config.yaml"
		" --force"
		" -o $output_dir"
		" -n $output_name"
		" {input}"
		" 1>{log} 2>&1"


