

rule mapping_DNA_pe:
	input:
		reads=rules.trimming_DNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
		idx_build=multiext("resources/{ref_name}/genome.fasta".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		), ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		bam=temp("results/{project}/mapping/dna/pe/{sample}-{unit}.sorted.bam"),
		idx=temp("results/{project}/mapping/dna/pe/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/{project}/log/mapping/dna/pe/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_pe"]["mapping_params"],
		sort_extra=config["mapping_DNA_pe"]["sort_params"],
		tmpdir=temp(directory("results/{project}/mapping/dna/pe/{sample}-{unit}.samtools_tmp")),
	threads: config["mapping_DNA_pe"]["threads"]
	conda:
		"../envs/bwa-mem2.yaml"
	shell:
		"("
		"bwa-mem2 mem"
		" -t {threads}"
		" {params.mapping_extra}"
		" {input.idx}"
		" {input.reads}"
		" |"
		" samtools sort"
		" -o {output.bam}"
		" --write-index"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
		" && rm -fr {params.tmpdir}"
		")"
		" 1>{log} 2>&1"


rule mapping_DNA_se:
	input:
		reads=rules.trimming_DNA_se.output.trimmed,
		idx="resources/{ref_name}/genome.fasta".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
		idx_build=multiext("resources/{ref_name}/genome.fasta".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		), ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		bam=temp("results/{project}/mapping/dna/se/{sample}-{unit}.sorted.bam"),
		idx=temp("results/{project}/mapping/dna/se/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/{project}/log/mapping/dna/se/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_se"]["mapping_params"],
		sort_extra=config["mapping_DNA_se"]["sort_params"],
		tmpdir=temp(directory("results/{project}/mapping/dna/pe/{sample}-{unit}.samtools_tmp")),
	threads: config["mapping_DNA_se"]["threads"]
	conda:
		"../envs/bwa-mem2.yaml"
	shell:
		"("
		"bwa-mem2 mem"
		" -t {threads}"
		" {params.mapping_extra}"
		" {input.idx}"
		" {input.reads}"
		" |"
		" samtools sort"
		" -o {output.bam}"
		" --write-index"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
		" && rm -fr {params.tmpdir}"
		")"
		" 1>{log} 2>&1"



rule mapping_RNA_pe:
	input:
		reads=rules.trimming_RNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta.STAR".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
	output:
		bam=temp("results/{project}/mapping/rna/pe/{sample}-{unit}.sorted.bam"),
		idx=temp("results/{project}/mapping/rna/pe/{sample}-{unit}.sorted.bam.csi"),
		tmpdir=temp(directory("results/{project}/mapping/rna/pe/{sample}-{unit}.STAR")),
	log:
		"results/{project}/log/mapping/rna/pe/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["mapping_RNA_pe"]["sjdbOverhang"],
		mapping_extra=config["mapping_RNA_pe"]["mapping_params"],
	threads: config["mapping_RNA_pe"]["threads"]
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
		" && mv {output.tmpdir}/STAR.Aligned.sortedByCoord.out.bam {output.bam}"
		" && samtools index -c {output.bam}"
		")"
		" 1>{log} 2>&1"


rule mapping_RNA_se:
	input:
		reads=rules.trimming_RNA_se.output.trimmed,
		idx="resources/{ref_name}/genome.fasta.STAR".format(
			ref_name=list(config["ref_genomes"].keys())[0],
		),
	output:
		bam=temp("results/{project}/mapping/rna/se/{sample}-{unit}.sorted.bam"),
		idx=temp("results/{project}/mapping/rna/se/{sample}-{unit}.sorted.bam.csi"),
		tmpdir=temp(directory("results/{project}/mapping/rna/se/{sample}-{unit}.STAR")),
	log:
		"results/{project}/log/mapping/rna/se/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["mapping_RNA_se"]["sjdbOverhang"],
		mapping_extra=config["mapping_RNA_se"]["mapping_params"],
	threads: config["mapping_RNA_se"]["threads"]
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
		" && mv {output.tmpdir}/STAR.Aligned.sortedByCoord.out.bam {output.bam}"
		" && samtools index -c {output.bam}"
		")"
		" 1>{log} 2>&1"


def get_bams_to_merge(wildcards):
	bams = []
	rows = samples.loc[(wildcards.sample), ["sample_id", "unit", "lib_type", "fq1", "fq2"]]
	for i, row in rows.iterrows():
		if pd.notnull(row.fq2):
			bams.append("results/{project}/mapping/{lib_type}/pe/{sample}-{unit}.sorted.bam".format(project=wildcards.project, sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		else:
			bams.append("results/{project}/mapping/{lib_type}/se/{sample}-{unit}.sorted.bam".format(project=wildcards.project, sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return bams


rule mapping_merge:
	input:
		unpack(get_bams_to_merge)
	output:
		bam="results/{project}/mapping_merged/{sample}.bam",
		idx="results/{project}/mapping_merged/{sample}.bam.csi",
	log:
		"results/{project}/log/mapping_merge/{sample}.log",
	params:
		extra=config["mapping_merge"]["params"],  # optional additional parameters, excluding --write-index which is implied by idx
	threads: config["mapping_merge"]["threads"]  # Samtools takes additional threads through its option -@
	conda:
		"../envs/samtools.yaml"
	shell:
		"("
		"samtools merge"
		" -o -"
		" {params.extra}"
		" {input} "
		" | samtools view "
		" -o {output.bam}"
		" --bam"
		" -F 4"
		" --write-index"
		")"
		" 1>{log} 2>&1"


