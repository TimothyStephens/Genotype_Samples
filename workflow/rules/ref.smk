

rule ref_parse_genome:
	output:
		"resources/{ref_name}/genome.fasta",
	log:
		"results/logs/resources/{ref_name}/ref_parse_genome.log",
	params:
		ref_file=lambda w: config["ref_genomes"][w.ref_name],
	conda:
		"../envs/bash.yaml"
	shell:
		"( if [[ {params.ref_file} == *.gz ]]; then zcat {params.ref_file} > {output}; else cat {params.ref_file} > {output}; fi ) 1>{log} 2>&1"


rule ref_faidx:
	input:
		"resources/{ref_name}/genome.fasta",
	output:
		"resources/{ref_name}/genome.fasta.fai",
	log:
		"results/logs/resources/{ref_name}/ref_faidx.log",
	conda:
		"../envs/bwa-mem2.yaml"
	shell:
		"samtools faidx"
		" {input}"
		" 1>{log} 2>&1"


rule ref_dict:
	input:
		"resources/{ref_name}/genome.fasta",
	output:
		"resources/{ref_name}/genome.dict",
	log:
		"results/logs/resources/{ref_name}/ref_dict.log",
	conda:
		"../envs/gatk4.yaml"
	shell:
		"("
		" export PATH=\"$CONDA_PREFIX/bin:$PATH\";"
		" gatk CreateSequenceDictionary -R {input} -O {output}"
		") 1>{log} 2>&1"


rule ref_DNA_mapping_index:
	input:
		"resources/{ref_name}/genome.fasta",
	output:
		multiext(
			"resources/{ref_name}/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"
		),
	log:
		"results/logs/resources/{ref_name}/ref_DNA_mapping_index.log",
	params:
		extra=config["ref_DNA_mapping_index"]["params"],
	conda:
		"../envs/bwa-mem2.yaml"
	shell:
		"bwa-mem2 index"
		" {input}"
		" {params.extra}"
		" 1>{log} 2>&1"


rule ref_RNA_mapping_index:
	input:
		"resources/{ref_name}/genome.fasta",
	output:
		directory("resources/{ref_name}/genome.fasta.STAR"),
	log:
		"results/logs/resources/{ref_name}/ref_RNA_mapping_index.log",
	params:
		extra=config["ref_RNA_mapping_index"]["params"],
		tmpdir=temp(directory("resources/{ref_name}/STARtmp")),
	threads: config["ref_RNA_mapping_index"]["threads"]
	conda:
		"../envs/star.yaml"
	shell:
		"("
		"rm -fr {params.tmpdir}; "
		"STAR"
		" --runThreadN {threads}"
		" --runMode genomeGenerate"
		" --genomeFastaFiles {input}"
		" {params.extra}"
		" --outTmpDir {params.tmpdir}"
		" --genomeDir {output}"
		" && rm -fr {params.tmpdir}"
		")"
		" 1>{log} 2>&1"


rule ref_DNA_longRead_mapping_index:
	input:
		"resources/{ref_name}/genome.fasta",
	output:
                "resources/{ref_name}/genome.fasta.hifi.mmi",
	log:
		"results/logs/resources/{ref_name}/ref_DNA_longRead_mapping_index.log",
	params:
		extra=config["ref_DNA_longRead_mapping_index"]["params"],
	conda:
		"../envs/pbmm2.yaml"
	shell:
		"pbmm2 index"
		" -j 1 --preset HiFi"
		" {params.extra}"
		" {input}"
		" {output}"
		" 1>{log} 2>&1"


rule ref_RNA_longRead_mapping_index:
	input:
		"resources/{ref_name}/genome.fasta",
	output:
		"resources/{ref_name}/genome.fasta.isoseq.mmi",
	log:
		"results/logs/resources/{ref_name}/ref_RNA_longRead_mapping_index.log",
	params:
		extra=config["ref_RNA_longRead_mapping_index"]["params"],
	conda:
		"../envs/pbmm2.yaml"
	shell:
		"pbmm2 index"
		" -j 1 --preset ISOSEQ"
		" {params.extra}"
		" {input}"
		" {output}"
		" 1>{log} 2>&1"

