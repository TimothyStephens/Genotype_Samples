

def get_fastq_DNA_pe(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "dna-pe"), ["fq1", "fq2"]]
	fq1 = fastqs.fq1.iloc[0]
	fq2 = fastqs.fq2.iloc[0]
	if fq1.startswith("DRR") or fq1.startswith("ERR") or fq1.startswith("SRR"):
		return {"sample": [
				"data/pe/{}_1.fastq.gz".format(fq1), 
				"data/pe/{}_2.fastq.gz".format(fq2),
				]}
	else:
		return {"sample": [fq1, fq2]}

def get_fastq_DNA_se(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "dna-se"), ["fq1"]]
	fq1 = fastqs.fq1.iloc[0]
	if fq1.startswith("DRR") or fq1.startswith("ERR") or fq1.startswith("SRR"):
		return {"sample": ["data/se/{}.fastq.gz".format(fq1)]}
	else:
		return {"sample": [fq1]}

def get_fastq_DNA_long(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "dna-long"), ["fq1"]]
	fq1 = fastqs.fq1.iloc[0]
	if fq1.startswith("DRR") or fq1.startswith("ERR") or fq1.startswith("SRR"):
		return {"sample": ["data/long/{}.fastq.gz".format(fq1)]}
	else:
		return {"sample": [fq1]}

def get_fastq_RNA_pe(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "rna-pe"), ["fq1", "fq2"]]
	fq1 = fastqs.fq1.iloc[0]
	fq2 = fastqs.fq2.iloc[0]
	if fq1.startswith("DRR") or fq1.startswith("ERR") or fq1.startswith("SRR"):
		return {"sample": [
				"data/pe/{}_1.fastq.gz".format(fq1), 
				"data/pe/{}_2.fastq.gz".format(fq2),
				]}
	else:
		return {"sample": [fq1, fq2]}

def get_fastq_RNA_se(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "rna-se"), ["fq1"]]
	fq1 = fastqs.fq1.iloc[0]
	if fq1.startswith("DRR") or fq1.startswith("ERR") or fq1.startswith("SRR"):
		return {"sample": ["data/se/{}.fastq.gz".format(fq1)]}
	else:
		return {"sample": [fq1]}

def get_fastq_RNA_long(wildcards):
	fastqs = samples.loc[(wildcards.sample, wildcards.unit, "rna-long"), ["fq1"]]
	fq1 = fastqs.fq1.iloc[0]
	if fq1.startswith("DRR") or fq1.startswith("ERR") or fq1.startswith("SRR"):
		return {"sample": ["data/long/{}.fastq.gz".format(fq1)]}
	else:
		return {"sample": [fq1]}


rule trimming_DNA_pe:
	input:
		unpack(get_fastq_DNA_pe),
	output:
		trimmed=[
			"results/trimmed/dna-pe/{sample}-{unit}.1.fastq.gz",
			"results/trimmed/dna-pe/{sample}-{unit}.2.fastq.gz",
		],
		html="results/qc/trimmed/dna-pe/{sample}-{unit}.html",
		json="results/qc/trimmed/dna-pe/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/dna-pe/{sample}-{unit}.log",
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
		trimmed=["results/trimmed/dna-se/{sample}-{unit}.1.fastq.gz"],
		html="results/qc/trimmed/dna-se/{sample}-{unit}.html",
		json="results/qc/trimmed/dna-se/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/dna-se/{sample}-{unit}.log",
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


rule trimming_DNA_long:
	input:
		unpack(get_fastq_DNA_long),
	output:
		trimmed=["results/trimmed/dna-long/{sample}-{unit}.1.fastq.gz"],
		html="results/qc/trimmed/dna-long/{sample}-{unit}.html",
		json="results/qc/trimmed/dna-long/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/dna-long/{sample}-{unit}.log",
	params:
		extra=config["trimming_DNA_long"]["params"],
	threads: config["trimming_DNA_long"]["threads"]
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
			"results/trimmed/rna-pe/{sample}-{unit}.1.fastq.gz",
			"results/trimmed/rna-pe/{sample}-{unit}.2.fastq.gz",
		],
		html="results/qc/trimmed/rna-pe/{sample}-{unit}.html",
		json="results/qc/trimmed/rna-pe/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/rna-pe/{sample}-{unit}.log",
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
		trimmed=["results/trimmed/rna-se/{sample}-{unit}.1.fastq.gz"],
		html="results/qc/trimmed/rna-se/{sample}-{unit}.html",
		json="results/qc/trimmed/rna-se/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/rna-se/{sample}-{unit}.log",
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


rule trimming_RNA_long:
	input:
		unpack(get_fastq_RNA_long),
	output:
		trimmed=["results/trimmed/rna-long/{sample}-{unit}.1.fastq.gz"],
		html="results/qc/trimmed/rna-long/{sample}-{unit}.html",
		json="results/qc/trimmed/rna-long/{sample}-{unit}_fastp.json",
	log:
		"results/logs/trimmed/rna-long/{sample}-{unit}.log",
	params:
                extra=config["trimming_RNA_long"]["params"],
	threads: config["trimming_RNA_long"]["threads"]
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


