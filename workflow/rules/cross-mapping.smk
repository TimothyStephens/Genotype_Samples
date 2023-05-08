

rule crossMapping_DNA_pe:
	input:
		reads=rules.trimming_DNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta",
		idx_build=multiext("resources/{ref_name}/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		"results/cross-mapping/{ref_name}/dna/pe/{sample}-{unit}.samtools_stats.txt"
	log:
		"results/logs/cross-mapping/{ref_name}/dna/pe/{sample}-{unit}.log",
	params:
		mapping_extra=config["crossMapping_DNA_pe"]["mapping_params"],
		sort_extra=config["crossMapping_DNA_pe"]["sort_params"],
		stats_extra=config["crossMapping_DNA_pe"]["stats_params"],
		tmpdir=temp(directory("results/cross-mapping/{ref_name}/dna/pe/{sample}-{unit}.samtools_tmp")),
	threads: config["crossMapping_DNA_pe"]["threads"]
	resources:
		mem_gb=config["crossMapping_DNA_pe"]["memory"]
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
		"results/cross-mapping/{ref_name}/dna/se/{sample}-{unit}.samtools_stats.txt"
	log:
		"results/logs/cross-mapping/{ref_name}/dna/se/{sample}-{unit}.log",
	params:
		mapping_extra=config["crossMapping_DNA_se"]["mapping_params"],
		sort_extra=config["crossMapping_DNA_se"]["sort_params"],
		stats_extra=config["crossMapping_DNA_se"]["stats_params"],
		tmpdir=temp(directory("results/cross-mapping/{ref_name}/dna/se/{sample}-{unit}.samtools_tmp")),
	threads: config["crossMapping_DNA_se"]["threads"]
	resources:
		mem_gb=config["crossMapping_DNA_se"]["memory"]
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
		stats="results/cross-mapping/{ref_name}/rna/pe/{sample}-{unit}.samtools_stats.txt",
		tmpdir=temp(directory("results/cross-mapping/{ref_name}/rna/pe/{sample}-{unit}.STAR")),
	log:
		"results/logs/cross-mapping/{ref_name}/rna/pe/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["crossMapping_RNA_pe"]["sjdbOverhang"],
		mapping_extra=config["crossMapping_RNA_pe"]["mapping_params"],
		stats_extra=config["crossMapping_RNA_pe"]["stats_params"],
	threads: config["crossMapping_RNA_pe"]["threads"]
	resources:
		mem_gb=config["crossMapping_RNA_pe"]["memory"]
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
		stats="results/cross-mapping/{ref_name}/rna/se/{sample}-{unit}.samtools_stats.txt",
		tmpdir=temp(directory("results/cross-mapping/{ref_name}/rna/se/{sample}-{unit}.STAR")),
	log:
		"results/logs/cross-mapping/{ref_name}/rna/se/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["crossMapping_RNA_se"]["sjdbOverhang"],
		mapping_extra=config["crossMapping_RNA_se"]["mapping_params"],
		stats_extra=config["crossMapping_RNA_se"]["stats_params"],
	threads: config["crossMapping_RNA_se"]["threads"]
	resources:
		mem_gb=config["crossMapping_RNA_se"]["memory"]
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
					ref_name=ref_name,
				))
			else:
				out.append("{ref_name}/{lib_type}/se/{sample}-{unit}".format(
						sample=row.sample_id, 
						unit=row.unit, 
						lib_type=row.lib_type,
						ref_name=ref_name,
					))
	return out


rule format_crossMapping_results:
	input:
		lambda wildcards: expand("results/cross-mapping/{s}.samtools_stats.txt", 
			s=expand_crossMapping_results_paths()),
	output:
		"results/{project}/final/mapping_rates.tsv"
	log:
		"results/logs/{project}/combine_mapping_results.log",
	script:
		"../scripts/format_crossMapping_results.sh"


