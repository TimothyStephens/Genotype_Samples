


def get_trimmed_fastqs(wildcards):
	fastqs = []
	rows = samples.loc[(wildcards.sample), ["sample_id", "unit", "lib_type", "fq1", "fq2"]]
	for i, row in rows.iterrows():
		if row.lib_type != 'dna':
			continue
		if pd.notnull(row.fq2):
			fastqs.append("results/trimmed/dna/pe/{sample}-{unit}.1.fastq.gz".format(sample=row.sample_id, unit=row.unit))
			fastqs.append("results/trimmed/dna/pe/{sample}-{unit}.2.fastq.gz".format(sample=row.sample_id, unit=row.unit))
		else:
			fastqs.append("results/trimmed/dna/se/{sample}-{unit}.1.fastq.gz".format(sample=row.sample_id, unit=row.unit))
	return(fastqs)


rule kmer_analysis_createDB_list_files:
	input:
		unpack(get_trimmed_fastqs),
	output:
		fof="results/kmer_analysis/kmerDB/{sample}.fof",
	log:
		"results/logs/kmer_analysis/kmerDB/{sample}.fof.log",
	conda:
		"../envs/bash.yaml"
	shell:
		"("
		"echo {input} | sed -e 's/ /\\n/g' > {output.fof}"
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_createDB:
	input:
		fof=rules.kmer_analysis_createDB_list_files.output.fof,
	output:
		db="results/kmer_analysis/kmerDB/{sample}.kcmdb",
		db_pre=temp("results/kmer_analysis/kmerDB/{sample}.kcmdb.kmc_pre"),
		db_suf=temp("results/kmer_analysis/kmerDB/{sample}.kcmdb.kmc_suf"),
		tmpdir=temp(directory("results/kmer_analysis/kmerDB/{sample}.kcmdb_tmp")),
	log:
		"results/logs/kmer_analysis/kmerDB/{sample}.create.log",
	params:
		kmer=config["kmer_analysis"]["kmer"],
		min_count=config["kmer_analysis"]["min_count"],
		max_count=config["kmer_analysis"]["max_count"],
		extra=config["kmer_analysis_createDB"]["params"],
	resources:
		mem_gb=config["kmer_analysis_createDB"]["memory"]
	threads: config["kmer_analysis_createDB"]["threads"]
	container:
		"docker://hamiltonjp/smudgeplot:0.2.4"
	shell:
		"("
		"mkdir -p {output.tmpdir}; "
		"kmc"
		" -k{params.kmer}"
		" -m{resources.mem_gb}"
		" -ci{params.min_count}"
		" -cs{params.max_count}"
		" -t{threads}"
		" -fq"
		" {params.extra}"
		" @{input.fof}"
		" {output.db}"
		" {output.tmpdir}"
		" && touch {output.db}"
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_histogram:
	input:
		db=rules.kmer_analysis_createDB.output.db,
		db_pre=rules.kmer_analysis_createDB.output.db_pre,
		db_suf=rules.kmer_analysis_createDB.output.db_suf,
	output:
		histo="results/kmer_analysis/histogram/{sample}.histo",
	log:
		"results/logs/kmer_analysis/histogram/{sample}.log",
	params:
		max_count=config["kmer_analysis"]["max_count"],
	container:
		"docker://hamiltonjp/smudgeplot:0.2.4"
	shell:
		"("
		"kmc_tools"
		" transform {input.db}"
		" histogram {output.histo} -cx{params.max_count}"
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_smudgeplot_cutoffs:
	input:
		histo=rules.kmer_analysis_histogram.output.histo,
	output:
		cutoffs="results/kmer_analysis/cutoffs/{sample}.kmer_cov",
		ploidy="results/kmer_analysis/cutoffs/{sample}.genomescope_ploidy",
	log:
		"results/logs/kmer_analysis/cutoffs/{sample}.cutoffs.log",
	container:
		"docker://hamiltonjp/smudgeplot:0.2.4"
	shell:
		"("
		"SMP=$(which smudgeplot.py); "
		"L=$(python $SMP cutoff {input.histo} L); "
		"U=$(python $SMP cutoff {input.histo} U); "
		"P=2; "
		"echo -e \"CUTOFF_LOWER\\t$L\" >  {output.cutoffs}; "
		"echo -e \"CUTOFF_UPPER\\t$U\" >> {output.cutoffs}; "
		"echo -e \"PLOIDY\\t$P\" >> {output.ploidy}; "
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_genomescope2:
	input:
		histo=rules.kmer_analysis_histogram.output.histo,
		ploidy=rules.kmer_analysis_smudgeplot_cutoffs.output.ploidy,
	output:
		plot_linear="results/{project}/kmer_analysis/genomescope2/{sample}_linear_plot.png",
		plot_log="results/{project}/kmer_analysis/genomescope2/{sample}_log_plot.png",
		plot_transformed_linear="results/{project}/kmer_analysis/genomescope2/{sample}_transformed_linear_plot.png",
		plot_transformed_log="results/{project}/kmer_analysis/genomescope2/{sample}_transformed_log_plot.png",
	log:
		"results/logs/{project}/kmer_analysis/genomescope2/{sample}.log",
	params:
		kmer=config["kmer_analysis"]["kmer"],
		extra=config["kmer_analysis_GenomeScope2"]["params"],
		outdir="results/{project}/kmer_analysis/genomescope2"
	container:
		"docker://abner12/genomescope:2.0"
	shell:
		"("
		"P=$(grep PLOIDY {input.ploidy} | cut -f2); "
		"genomescope.R"
		" --ploidy $P"
		" -i {input.histo}"
		" -o {params.outdir}"
		" -n {wildcards.sample}"
		" -k {params.kmer}"
		" {params.extra}"
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_smudgeplot_tranform:
	input:
		cutoffs=rules.kmer_analysis_smudgeplot_cutoffs.output.cutoffs,
		db=rules.kmer_analysis_createDB.output.db,
		db_pre=rules.kmer_analysis_createDB.output.db_pre,
		db_suf=rules.kmer_analysis_createDB.output.db_suf,
	output:
		db_reduced="results/kmer_analysis/kmerDB/{sample}.reduce.kcmdb",
		db_reduced_pre=temp("results/kmer_analysis/kmerDB/{sample}.reduce.kcmdb.kmc_pre"),
		db_reduced_suf=temp("results/kmer_analysis/kmerDB/{sample}.reduce.kcmdb.kmc_suf"),
	log:
		"results/logs/kmer_analysis/kmerDB/{sample}.reduce.log",
	container:
		"docker://hamiltonjp/smudgeplot:0.2.4"
	shell:
		"("
		"L=$(grep CUTOFF_LOWER {input.cutoffs} | cut -f2); "
		"U=$(grep CUTOFF_UPPER {input.cutoffs} | cut -f2); "
		"kmc_tools"
		" transform {input.db} -ci$L -cx$U"
		" reduce {output.db_reduced}"
		" && touch {output.db_reduced}"
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_smudgeplot_pairs:
	input:
		db_reduced=rules.kmer_analysis_smudgeplot_tranform.output.db_reduced,
		db_reduced_pre=rules.kmer_analysis_smudgeplot_tranform.output.db_reduced_pre,
		db_reduced_suf=rules.kmer_analysis_smudgeplot_tranform.output.db_reduced_suf,
	output:
		cov="results/kmer_analysis/smudgeplot/{sample}.smudge_pairs.coverages.tsv",
		pairs="results/kmer_analysis/smudgeplot/{sample}.smudge_pairs.pairs.tsv",
		family="results/kmer_analysis/smudgeplot/{sample}.smudge_pairs.familysizes.tsv",
	log:
		cov="results/logs/kmer_analysis/smudgeplot/{sample}.smudge_pairs.log",
	resources:
		mem_gb=config["kmer_analysis_smudgeplot_pairs"]["memory"]
	container:
		"docker://hamiltonjp/smudgeplot:0.2.4"
	shell:
		"("
		"smudge_pairs"
		" {input.db_reduced}"
		" {output.cov}"
		" {output.pairs}"
		" > {output.family}"
		")"
		" 1>{log} 2>&1"


rule kmer_analysis_smudgeplot:
	input:
		cov=rules.kmer_analysis_smudgeplot_pairs.output.cov,
	output:
		plot="results/{project}/kmer_analysis/smudgeplot/{sample}_smudgeplot.png",
		plot_log="results/{project}/kmer_analysis/smudgeplot/{sample}_smudgeplot_log10.png",
	params:
		kmer=config["kmer_analysis"]["kmer"],
		extra=config["kmer_analysis_Smudgeplot"]["params"],
		out_prefix="results/{project}/kmer_analysis/smudgeplot/{sample}",
	log:
		"results/logs/{project}/kmer_analysis/smudgeplot/{sample}_smudgeplot.log",
	container:
		"docker://hamiltonjp/smudgeplot:0.2.4"
	shell:
		"("
		"SMP=$(which smudgeplot.py); "
		"python $SMP plot"
		" -k {params.kmer}"
		" -o {params.out_prefix}"
		" {params.extra}"
		" {input.cov}"
		")"
		" 1>{log} 2>&1"


