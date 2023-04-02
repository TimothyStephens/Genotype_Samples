

def get_fastq_DNA_pe(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "dna"), ["fq1", "fq2"]]
	if fastqs.fq1.startswith("SRR"):
		return {"sample": ["data/pe/"+fastqs.fq1+"_1.fastq.gz", "data/pe/"+fastqs.fq2+"_2.fastq.gz"]}
	else:
		return {"sample": [fastqs.fq1, fastqs.fq2]}

def get_fastq_DNA_se(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "dna"), ["fq1"]]
	if fastqs.fq1.startswith("SRR"):
		return {"sample": ["data/se/"+fastqs.fq1+".fastq.gz"]}
	else:
		return {"sample": [fastqs.fq1]}

def get_fastq_RNA_pe(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "rna"), ["fq1", "fq2"]]
	if fastqs.fq1.startswith("SRR"):
		return {"sample": ["data/pe/"+fastqs.fq1+"_1.fastq.gz", "data/pe/"+fastqs.fq2+"_2.fastq.gz"]}
	else:
		return {"sample": [fastqs.fq1, fastqs.fq2]}

def get_fastq_RNA_se(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "rna"), ["fq1"]]
	if fastqs.fq1.startswith("SRR"):
		return {"sample": ["data/se/"+fastqs.fq1+".fastq.gz"]}
	else:
		return {"sample": [fastqs.fq1]}


rule trimming_DNA_pe:
	input:
		unpack(get_fastq_DNA_pe),
	output:
		trimmed=[
			"results/trimmed/dna/pe/{sample}-{unit}.1.fastq.gz",
			"results/trimmed/dna/pe/{sample}-{unit}.2.fastq.gz",
		],
		html="results/qc/trimmed/dna/pe/{sample}-{unit}.html",
		json="results/qc/trimmed/dna/pe/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/dna/pe/{sample}-{unit}.log",
	params:
		extra=config["trimming_DNA_pe"]["params"],
	threads: config["trimming_DNA_pe"]["threads"]
	conda:
		"../envs/fastp.yaml"
	shell:
		"fastp"
		" --in1 {input[0]}"
		" --in2 {input[1]}"
		" --out1 {output.trimmed[0]}"
		" --out2 {output.trimmed[1]}"
		" --json {output.json}"
		" --html {output.html}"
		" --thread {threads}"
		" {params.extra}"
		" 1>{log} 2>&1"



rule trimming_DNA_se:
	input:
		unpack(get_fastq_DNA_se),
	output:
		trimmed=["results/trimmed/dna/se/{sample}-{unit}.1.fastq.gz"],
		html="results/qc/trimmed/dna/se/{sample}-{unit}.html",
		json="results/qc/trimmed/dna/se/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimming/dna/se/{sample}-{unit}.log",
	params:
		extra=config["trimming_DNA_se"]["params"],
	threads: config["trimming_DNA_se"]["threads"]
	conda:
		"../envs/fastp.yaml"
	shell:
		"fastp"
		" --in1 {input[0]}"
		" --out1 {output.trimmed[0]}"
		" --json {output.json}"
		" --html {output.html}"
		" --thread {threads}"
		" {params.extra}"
		" 1>{log} 2>&1"


rule trimming_RNA_pe:
	input:
		unpack(get_fastq_RNA_pe),
	output:
		trimmed=[
			"results/trimmed/rna/pe/{sample}-{unit}.1.fastq.gz",
			"results/trimmed/rna/pe/{sample}-{unit}.2.fastq.gz",
		],
		html="results/qc/trimmed/rna/pe/{sample}-{unit}.html",
		json="results/qc/trimmed/rna/pe/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/rna/pe/{sample}-{unit}.log",
	params:
		extra=config["trimming_RNA_pe"]["params"],
	threads: config["trimming_RNA_pe"]["threads"]
	conda:
		"../envs/fastp.yaml"
	shell:
		"fastp"
		" --in1 {input[0]}"
		" --in2 {input[1]}"
		" --out1 {output.trimmed[0]}"
		" --out2 {output.trimmed[1]}"
		" --json {output.json}"
		" --html {output.html}"
		" --thread {threads}"
		" {params.extra}"
		" 1>{log} 2>&1"


rule trimming_RNA_se:
	input:
		unpack(get_fastq_RNA_se),
	output:
		trimmed=["results/trimmed/rna/se/{sample}-{unit}.1.fastq.gz"],
		html="results/qc/trimmed/rna/se/{sample}-{unit}.html",
		json="results/qc/trimmed/rna/se/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimming/rna/se/{sample}-{unit}.log",
	params:
		extra=config["trimming_RNA_se"]["params"],
	threads: config["trimming_RNA_se"]["threads"]
	conda:
		"../envs/fastp.yaml"
	shell:
		"fastp"
		" --in1 {input[0]}"
		" --out1 {output.trimmed[0]}"
		" --json {output.json}"
		" --html {output.html}"
		" --thread {threads}"
		" {params.extra}"
		" 1>{log} 2>&1"


