

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
		if wildcards.fq == "{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit):
			if row.fq1.startswith("SRR") and pd.notnull(row.fq2):
				return ["data/pe/{}_1.fastq.gz".format(row.fq1)]
			elif row.fq1.startswith("SRR") and pd.isnull(row.fq2):
				return ["data/se/{}.fastq.gz".format(row.fq1)]
			else:
				return [row.fq1]
		elif wildcards.fq == "{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit):
			if row.fq2.startswith("SRR"):
				return ["data/pe/{}_2.fastq.gz".format(row.fq2)]
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
		bam="results/mapping_merged/{ref_name}/{sample}.bam",
	output:
		"results/qc/mapping_merged/{ref_name}/{sample}_samtools_stats.txt",
	log:
		"results/logs/qc/mapping_merged/{ref_name}/{sample}.log",
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
		bam="results/mapping_merged/{ref_name}/{sample}.bam",
	output:
		directory("results/qc/mapping_merged/{ref_name}/{sample}_qualimap"),
	log:
		"results/logs/qc/mapping_merged/{ref_name}/{sample}.qualimap.log",
	params:
		extra=config["qc_qualimap"]["params"],
	resources:
		mem_gb=config["qc_qualimap"]["memory"]
	conda:
		"../envs/qualimap.yaml"
	shell:
		"qualimap bamqc"
		" {params.extra}"
		" --java-mem-size={resources.mem_gb}G"
		" -bam {input.bam}"
		" -outdir {output}"
		" -outformat HTML"
		" 1>{log} 2>&1"


rule qc_samtools_stats_unmerged:
	input:
		bam="results/mapping/{ref_name}/{fq}.sorted.bam",
	output:
		"results/qc/mapping/{ref_name}/{fq}_samtools_stats.txt",
	log:
		"results/logs/qc/mapping/{ref_name}/{fq}.log",
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
		"results/{project}/calling_merged/{name}.vcf.gz",
	output:
		"results/{project}/qc/calling_merged/{name}_bcftools_stats.txt",
	log:
		"results/logs/{project}/qc/calling_merged/{name}_bcftools_stats.log",
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


