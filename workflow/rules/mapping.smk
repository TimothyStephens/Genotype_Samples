

rule mapping_DNA_pe:
	input:
		reads=rules.trimming_DNA_pe.output.trimmed,
		idx="resources/{ref_name}/genome.fasta",
		idx_build=multiext("resources/{ref_name}/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
	output:
		bam=temp("results/mapping/{ref_name}/dna-pe/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapping/{ref_name}/dna-pe/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/logs/mapping/{ref_name}/dna-pe/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_pe"]["mapping_params"],
		sort_extra=config["mapping_DNA_pe"]["sort_params"],
		tmpdir=temp(directory("results/mapping/{ref_name}/dna-pe/{sample}-{unit}.samtools_tmp")),
	threads: config["mapping_DNA_pe"]["threads"]
	resources:
		mem_gb=config["mapping_DNA_pe"]["memory"]
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
		bam=temp("results/mapping/{ref_name}/dna-se/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapping/{ref_name}/dna-se/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/logs/mapping/{ref_name}/dna-se/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_se"]["mapping_params"],
		sort_extra=config["mapping_DNA_se"]["sort_params"],
		tmpdir=temp(directory("results/mapping/{ref_name}/dna-pe/{sample}-{unit}.samtools_tmp")),
	threads: config["mapping_DNA_se"]["threads"]
	resources:
		mem_gb=config["mapping_DNA_se"]["memory"]
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


rule mapping_DNA_long:
	input:
		reads=rules.trimming_DNA_long.output.trimmed,
		idx="resources/{ref_name}/genome.fasta",
		idx_build="resources/{ref_name}/genome.fasta.hifi.mmi",
	output:
		fofn="results/mapping/{ref_name}/dna-long/{sample}-{unit}.fofn",
		bam=temp("results/mapping/{ref_name}/dna-long/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapping/{ref_name}/dna-long/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/logs/mapping/{ref_name}/dna-long/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_DNA_long"]["mapping_params"],
		sort_extra=config["mapping_DNA_long"]["sort_params"],
		tmpdir=temp(directory("results/mapping/{ref_name}/dna-long/{sample}-{unit}.samtools_tmp")),
	threads: config["mapping_DNA_long"]["threads"]
	resources:
		mem_gb=config["mapping_DNA_long"]["memory"]
	conda:
		"../envs/pbmm2.yaml"
	shell:
		"("
		"echo {input.reads} | sed -e 's/ /\\n/g' > {output.fofn}; "
		"pbmm2 align"
		" --preset HiFi"
		" --unmapped"
		" -j {threads}"
		" {params.mapping_extra}"
		" {input.idx_build}"
		" {output.fofn}"
		" | samtools sort"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
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
		bam=temp("results/mapping/{ref_name}/rna-pe/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapping/{ref_name}/rna-pe/{sample}-{unit}.sorted.bam.csi"),
		tmpdir=temp(directory("results/mapping/{ref_name}/rna-pe/{sample}-{unit}.STAR")),
	log:
		"results/logs/mapping/{ref_name}/rna-pe/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["mapping_RNA_pe"]["sjdbOverhang"],
		mapping_extra=config["mapping_RNA_pe"]["mapping_params"],
	threads: config["mapping_RNA_pe"]["threads"]
	resources:
		mem_gb=config["mapping_RNA_pe"]["memory"]
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
		bam=temp("results/mapping/{ref_name}/rna-se/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapping/{ref_name}/rna-se/{sample}-{unit}.sorted.bam.csi"),
		tmpdir=temp(directory("results/mapping/{ref_name}/rna-se/{sample}-{unit}.STAR")),
	log:
		"results/logs/mapping/{ref_name}/rna-se/{sample}-{unit}.log",
	params:
		sjdbOverhang=config["mapping_RNA_se"]["sjdbOverhang"],
		mapping_extra=config["mapping_RNA_se"]["mapping_params"],
	threads: config["mapping_RNA_se"]["threads"]
	resources:
		mem_gb=config["mapping_RNA_se"]["memory"]
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


rule mapping_RNA_long:
	input:
		reads=rules.trimming_RNA_long.output.trimmed,
		idx="resources/{ref_name}/genome.fasta",
		idx_build="resources/{ref_name}/genome.fasta.isoseq.mmi",
	output:
		fofn="results/mapping/{ref_name}/rna-long/{sample}-{unit}.fofn",
		bam=temp("results/mapping/{ref_name}/rna-long/{sample}-{unit}.sorted.bam"),
		idx=temp("results/mapping/{ref_name}/rna-long/{sample}-{unit}.sorted.bam.csi"),
	log:
		"results/logs/mapping/{ref_name}/rna-long/{sample}-{unit}.log",
	params:
		mapping_extra=config["mapping_RNA_long"]["mapping_params"],
		sort_extra=config["mapping_RNA_long"]["sort_params"],
		tmpdir=temp(directory("results/mapping/{ref_name}/rna-long/{sample}-{unit}.samtools_tmp")),
	threads: config["mapping_RNA_long"]["threads"]
	resources:
		mem_gb=config["mapping_RNA_long"]["memory"]
	conda:
		"../envs/pbmm2.yaml"
	shell:
		"("
		"echo {input.reads} | sed -e 's/ /\\n/g' > {output.fofn}; "
		"pbmm2 align"
		" --preset ISOSEQ"
		" --unmapped"
		" -j {threads}"
		" {params.mapping_extra}"
		" {input.idx_build}"
		" {output.fofn}"
		" | samtools sort"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
		" |"
		" samtools sort"
		" -o {output.bam}"
		" --write-index"
		" {params.sort_extra}"
		" -T {params.tmpdir}"
		" && rm -fr {params.tmpdir}"
		")"
		" 1>{log} 2>&1"


def get_bams_to_merge(wildcards):
	bams = []
	rows = samples.loc[(wildcards.sample), ["sample_id", "unit", "lib_type", "fq1", "fq2"]]
	for i, row in rows.iterrows():
		bams.append("results/mapping/{ref_name}/{lib_type}/{sample}-{unit}.sorted.bam".format(ref_name=wildcards.ref_name, sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return bams


rule mapping_merge:
	input:
		unpack(get_bams_to_merge)
	output:
		bam="results/mapping_merged/{ref_name}/{sample}.bam",
		idx="results/mapping_merged/{ref_name}/{sample}.bam.csi",
	log:
		"results/logs/mapping_merge/{ref_name}/{sample}.log",
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


