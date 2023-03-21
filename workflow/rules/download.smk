

rule download_fastq_pe:
	output:
		temp("data/pe/{accession}_1.fastq.gz"),
		temp("data/pe/{accession}_2.fastq.gz"),
	log:
		"logs/download_fastq_pe/{accession}.gz.log"
	params:
		extra="--skip-technical " + config["download_fastq"]["params"]
	threads: config["download_fastq"]["threads"]
	wrapper:
		"v1.23.4/bio/sra-tools/fasterq-dump"


rule download_fastq_se:
	output:
		temp("data/se/{accession}.fastq.gz")
	log:
		"logs/download_fastq_se/{accession}.gz.log"
	params:
		extra="--skip-technical " + config["download_fastq"]["params"]
	threads: config["download_fastq"]["threads"]
	wrapper:
		"v1.23.4/bio/sra-tools/fasterq-dump"


