

rule relatedness_vcftools_relatedness2:
	input:
		rules.calling_filter_merged_VCF.output,
	output:
		"results/{project}/relatedness/calls.filtered.vcf.relatedness2",
	log:
		"results/logs/{project}/relatedness/vcftools_relatedness2.log",
	params:
		out_prefix="results/{project}/relatedness/calls.filtered.vcf",
		extra=config["relatedness_vcftools_relatedness2"]["params"],
	conda:
		"../envs/vcftools.yaml"
	shell:
		"vcftools "
		"--gzvcf {input} "
		"--out {params.out_prefix} "
		"{params.extra} "
		"--relatedness2 "
		"1>{log} 2>&1"


rule relatedness_vcf_clone_detect_count:
	input:
		rules.calling_filter_merged_VCF.output,
	output:
		"results/{project}/relatedness/calls.filtered.vcf.allelic_similarity.csv"
	log:
		"results/logs/{project}/relatedness/vcf_clone_detect-count.log"
	conda:
		"../envs/vcf_clone_detect.yaml"
	shell:
		"python workflow/scripts/vcf_clone_detect-count.py"
		" --vcf <(gunzip -c {input})"
		" --output {output}"
		" 1>{log} 2>&1"


rule relatedness_vcf_clone_detect_group:
	input:
		rules.relatedness_vcf_clone_detect_count.output,
	output:
		"results/{project}/relatedness/calls.filtered.vcf.allelic_similarity.groups.csv"
	log:
		"results/logs/{project}/relatedness/vcf_clone_detect-groups.log"
	params:
		threshold=config["relatedness_vcf_clone_detect"]["threshold"],
	conda:
		"../envs/vcf_clone_detect.yaml"
	shell:
		"python workflow/scripts/vcf_clone_detect-group.py"
		" --input {input}"
		" --output {output}"
		" --threshold {params.threshold}"
		" 1>{log} 2>&1"


rule relatedness_write_bam_list:
	input:
		lambda wildcards: expand("results/mapping_merged/{ref_name}/{sample}.bam", 
			ref_name=list(config["ref_genomes"].keys())[0],
			sample=samples.sample_id.unique())
	output:
		files="results/{project}/relatedness/bam.filelist",
		labels="results/{project}/relatedness/bam.filelist.labels",
	log:
		"results/logs/{project}/relatedness/write_bam_list.log",
	conda:
		"../envs/bash.yaml"
	shell:'''
		rm -fr {output.files} {output.labels}
		for f in {input};
		do
		  echo $f >> {output.files}
		  echo $f | sed -e 's@.*/@@' -e 's/.bam//' >> {output.labels}
		done
		'''


rule relatedness_ANGSD_for_NgsRelate:
	input:
		rules.relatedness_write_bam_list.output.files,
	output:
		mafs="results/{project}/relatedness/NgsRelate.angsd.mafs.gz",
		glf="results/{project}/relatedness/NgsRelate.angsd.glf.gz",
	log:
		"results/logs/{project}/relatedness/NgsRelate.angsd.log",
	params:
		out_prefix="results/{project}/relatedness/NgsRelate.angsd",
		extra=config["relatedness_ANGSD_for_NgsRelate"]["params"],
	threads: config["relatedness_ANGSD_for_NgsRelate"]["threads"]
	container:
		"docker://lifebitai/angsd:0.933"
	shell:
		"(angsd -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -nThreads {threads} {params.extra} -bam {input} -out {params.out_prefix}) 1>{log} 2>&1"


rule relatedness_NgsRelate:
	input:
		mafs=rules.relatedness_ANGSD_for_NgsRelate.output.mafs,
		glf=rules.relatedness_ANGSD_for_NgsRelate.output.glf,
		labels=rules.relatedness_write_bam_list.output.labels,
	output:
		"results/{project}/relatedness/NgsRelate.results.tsv",
	log:
		"results/logs/{project}/relatedness/NgsRelate.log",
	params:
		extra=config["relatedness_NgsRelate"]["params"],
	threads: config["relatedness_NgsRelate"]["threads"]
	container:
		"docker://didillysquat/ngsrelate:latest"
	shell:
		"(zcat {input.mafs} | cut -f5 |sed 1d > {input.mafs}.freq; "
		"N=$(cat {input.labels} | wc -l); "
		"ngsRelate -g {input.glf} -n $N -f {input.mafs}.freq -z {input.labels} -O {output} -p {threads} {params.extra}) 1>{log} 2>&1"


rule relatedness_ANGSD_for_PCAngsd:
	input:
		rules.relatedness_write_bam_list.output.files,
	output:
		mafs="results/{project}/relatedness/PCAngsd.angsd.mafs.gz",
		beagle="results/{project}/relatedness/PCAngsd.angsd.beagle.gz",
	log:
		"results/logs/{project}/relatedness/PCAngsd.angsd.log",
	params:
		out_prefix="results/{project}/relatedness/PCAngsd.angsd",
		extra=config["relatedness_ANGSD_for_PCAngsd"]["params"],
	threads: config["relatedness_ANGSD_for_PCAngsd"]["threads"]
	container:
		"docker://lifebitai/angsd:0.933"
	shell:
		"(angsd -GL 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -nThreads {threads} {params.extra} -bam {input} -out {params.out_prefix}) 1>{log} 2>&1"


rule relatedness_PCAngsd_IndAlleleFreq:
	input:
		beagle=rules.relatedness_ANGSD_for_PCAngsd.output.beagle,
	output:
		"results/{project}/relatedness/PCAngsd.IndAlleleFreq.cov",
	log:
		"results/logs/{project}/relatedness/PCAngsd.IndAlleleFreq.log",
	params:
		out_prefix="results/{project}/relatedness/PCAngsd.IndAlleleFreq",
		extra=config["relatedness_PCAngsd_IndAlleleFreq"]["params"],
	threads: config["relatedness_PCAngsd_IndAlleleFreq"]["threads"]
	container:
		"docker://zjnolen/pcangsd:latest"
	shell:
		"(pcangsd --threads {threads} {params.extra} --beagle {input.beagle} --out {params.out_prefix}) 1>{log} 2>&1"


rule relatedness_PCAngsd_WithOutIndAlleleFreq:
	input:
		beagle=rules.relatedness_ANGSD_for_PCAngsd.output.beagle,
	output:
		"results/{project}/relatedness/PCAngsd.WithOutIndAlleleFreq.cov",
	log:
		"results/logs/{project}/relatedness/PCAngsd.WithOutIndAlleleFreq.log",
	params:
		out_prefix="results/{project}/relatedness/PCAngsd.WithOutIndAlleleFreq",
		extra=config["relatedness_PCAngsd_WithOutIndAlleleFreq"]["params"],
	threads: config["relatedness_PCAngsd_WithOutIndAlleleFreq"]["threads"]
	container:
		"docker://zjnolen/pcangsd:latest"
	shell:
		"(pcangsd --threads {threads} {params.extra} --beagle {input.beagle} --out {params.out_prefix} --iter 0) 1>{log} 2>&1"


rule relatedness_PCAngsd_Admixture:
	input:
		beagle=rules.relatedness_ANGSD_for_PCAngsd.output.beagle,
	output:
		"results/{project}/relatedness/PCAngsd.Admixture.admix.Q",
	log:
		"results/logs/{project}/relatedness/PCAngsd.Admixture.log",
	params:
		out_prefix="results/{project}/relatedness/PCAngsd.Admixture",
		extra=config["relatedness_PCAngsd_Admixture"]["params"],
	threads: config["relatedness_PCAngsd_Admixture"]["threads"]
	container:
		"docker://zjnolen/pcangsd:latest"
	shell:
		"(pcangsd --threads {threads} {params.extra} --beagle {input.beagle} --out {params.out_prefix} --admix --admix_alpha 50; cp {params.out_prefix}.admix.*.Q {params.out_prefix}.admix.Q) 1>{log} 2>&1"


rule format_results_PCAngsd_IndAlleleFreq:
	input:
		rules.relatedness_write_bam_list.output.labels,
		rules.relatedness_PCAngsd_IndAlleleFreq.output,
	output:
		"results/{project}/final/PCAngsd.IndAlleleFreq.tsv",
	log:
		"results/logs/{project}/relatedness/format_PCAngsd_IndAlleleFreq_results.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_PCAngsd_results.sh"


rule format_results_PCAngsd_WithOutIndAlleleFreq:
	input:
		rules.relatedness_write_bam_list.output.labels,
		rules.relatedness_PCAngsd_WithOutIndAlleleFreq.output,
	output:
		"results/{project}/final/PCAngsd.WithOutIndAlleleFreq.tsv",
	log:
		"results/logs/{project}/relatedness/format_PCAngsd_WithOutIndAlleleFreq_results.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_PCAngsd_results.sh"


rule format_results_PCAngsd_Admixture:
	input:
		rules.relatedness_write_bam_list.output.labels,
		rules.relatedness_PCAngsd_Admixture.output,
	output:
		"results/{project}/final/PCAngsd.Admixture.tsv",
	log:
		"results/logs/{project}/relatedness/format_PCAngsd_Admixture_results.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_PCAngsd_Admixture_results.sh"


rule format_results_vcftools_relatedness2:
	input:
		rules.relatedness_write_bam_list.output.labels,
		rules.relatedness_vcftools_relatedness2.output,
	output:
		"results/{project}/final/vcftools_relatedness2.tsv",
	log:
		"results/logs/{project}/relatedness/format_vcftools_relatedness2_results.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_vcftools_relatedness2_results.sh"


rule format_results_vcf_clone_detect_counts:
	input:
		rules.relatedness_write_bam_list.output.labels,
		rules.relatedness_vcf_clone_detect_count.output,
	output:
		"results/{project}/final/vcf_clone_detect.counts.tsv",
	log:
		"results/logs/{project}/relatedness/format_vcf_clone_detect_counts.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_vcf_clone_detect_counts.sh"


rule format_results_vcf_clone_detect_groups:
	input:
		rules.relatedness_vcf_clone_detect_group.output,
	output:
		"results/{project}/final/vcf_clone_detect.groups.tsv",
	log:
		"results/logs/{project}/relatedness/format_vcf_clone_detect_groups.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_vcf_clone_detect_groups.sh"


rule format_results_NgsRelate:
	input:
		rules.relatedness_NgsRelate.output,
	output:
		"results/{project}/final/NgsRelate.tsv",
	log:
		"results/logs/{project}/relatedness/format_NgsRelate_results.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_NgsRelate_results.sh"


