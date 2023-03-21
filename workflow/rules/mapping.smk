

rule mapping_DNA_pe:
	input:
		reads=rules.trimming_DNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta".format(ref_name=config["ref"]["name"]),
		idx_build=multiext("resources/{ref_name}/genome.fasta".format(ref_name=config["ref"]["name"]), ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		bam=temp("results/mapped/dna/pe/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapped/dna/pe/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/logs/mapped/dna/pe/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_pe"]["mapping_params"],
		sort_extra=config["mapping_DNA_pe"]["sort_params"],
		tmpdir=temp(directory("results/logs/mapped/dna/pe/{sample}-{unit}.samtools_tmp")),
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
		idx="resources/{ref_name}/genome.fasta".format(ref_name=config["ref"]["name"]),
		idx_build=multiext("resources/{ref_name}/genome.fasta".format(ref_name=config["ref"]["name"]), ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		bam=temp("results/mapped/dna/se/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapped/dna/se/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/logs/mapped/dna/se/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_se"]["mapping_params"],
		sort_extra=config["mapping_DNA_se"]["sort_params"],
		tmpdir=temp(directory("results/logs/mapped/dna/pe/{sample}-{unit}.samtools_tmp")),
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
		idx="resources/{ref_name}/genome.fasta.STAR".format(ref_name=config["ref"]["name"]),
	output:
		bam=temp("results/mapped/rna/pe/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapped/rna/pe/{sample}-{unit}.sorted.bam.csi"),
		tmpdir=temp(directory("results/mapped/rna/pe/{sample}-{unit}.STAR")),
	log:
		"results/logs/mapped/rna/pe/{sample}-{unit}.log",
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
		" --twopassMode Basic"
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
		idx="resources/{ref_name}/genome.fasta.STAR".format(ref_name=config["ref"]["name"]),
	output:
		bam=temp("results/mapped/rna/se/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapped/rna/se/{sample}-{unit}.sorted.bam.csi"),
		tmpdir=temp(directory("results/mapped/rna/se/{sample}-{unit}.STAR")),
	log:
		"results/logs/mapped/rna/se/{sample}-{unit}.log",
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
		" --twopassMode Basic"
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
		if row.lib_type == "dna" and pd.notnull(row.fq2):
			bams.append("results/mapped/{lib_type}/pe/{sample}-{unit}.sorted.bam".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "dna" and pd.isnull(row.fq2):
			bams.append("results/mapped/{lib_type}/se/{sample}-{unit}.sorted.bam".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "rna" and pd.notnull(row.fq2):
			bams.append("results/mapped/{lib_type}/pe/{sample}-{unit}.sorted.bam".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if row.lib_type == "rna" and pd.isnull(row.fq2):
			bams.append("results/mapped/{lib_type}/se/{sample}-{unit}.sorted.bam".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return list(set(bams))


rule mapping_merge:
	input:
		unpack(get_bams_to_merge)
	output:
		bam="results/mapped/{sample}.bam",
		idx="results/mapped/{sample}.bam.csi",
	log:
		"results/logs/samtools_merge/{sample}.log",
	params:
		extra=config["mapping_merge"]["params"],  # optional additional parameters, excluding --write-index which is implied by idx
	threads: config["mapping_merge"]["threads"]  # Samtools takes additional threads through its option -@
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools merge"
		" -o {output.bam}"
		" --write-index"
		" {params.extra}"
		" {input}"
		" 1>{log} 2>&1"


