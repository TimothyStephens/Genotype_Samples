

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
	resources:
		mem_gb=config["relatedness_ANGSD_for_NgsRelate"]["memory"]
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


rule relatedness_plink_PCA:
	input:
		vcf=rules.calling_filter_merged_VCF.output,
	output:
		eigenval="results/{project}/relatedness/plink.PCA.eigenval",
		eigenvec="results/{project}/relatedness/plink.PCA.eigenvec",
		nosex="results/{project}/relatedness/plink.PCA.nosex",
	log:
		"results/logs/{project}/relatedness/plink.PCA.log",
	params:
		out="results/{project}/relatedness/plink.PCA"
	conda:
		"../envs/plink.yaml"
	shell:
		"("
		"plink --vcf {input.vcf} --double-id --vcf-half-call m --allow-extra-chr"
		" --set-missing-var-ids @:#"
		" --indep-pairwise 50 10 0.1"
		" --out {params.out}; "
		"plink --vcf {input.vcf} --double-id --vcf-half-call m --allow-extra-chr"
		" --set-missing-var-ids @:#"
		" --extract {params.out}.prune.in"
		" --make-bed --pca --out {params.out}"
		") 1>{log} 2>&1"


rule relatedness_plink_Admixture_prepData:
	input:
		vcf=rules.calling_filter_merged_VCF.output,
	output:
		"results/{project}/relatedness/plink.Admixture.bed",
	log:
		"results/logs/{project}/relatedness/plink.Admixture.prepData.log",
	params:
		out="results/{project}/relatedness/plink.Admixture",
	conda:
		"../envs/plink.yaml"
	shell:
		"("
		"plink --vcf {input.vcf} --double-id --vcf-half-call m --allow-extra-chr --make-bed --out {params.out}; "
		"awk '{{$1=\"0\"; print $0}}' {params.out}.bim > {params.out}.bim.tmp; "
		"mv {params.out}.bim.tmp {params.out}.bim"
		") 1>{log} 2>&1"


rule relatedness_plink_Admixture:
	input:
		rules.relatedness_plink_Admixture_prepData.output,
	output:
		best="results/{project}/relatedness/plink.Admixture.best.Q",
	log:
		"results/logs/{project}/relatedness/plink.Admixture.log",
	params:
		out_dir="results/{project}/relatedness",
		out_prefix="plink.Admixture",
		Kmin=config["relatedness_plink_Admixture"]["Kmin"],
		Kmax=config["relatedness_plink_Admixture"]["Kmax"],
	conda:
		"../envs/admixture.yaml"
	threads: config["relatedness_plink_Admixture"]["threads"]
	shell:
		"("
		" ( cd {params.out_dir}; for i in {{{params.Kmin}..{params.Kmax}}}; do echo \"admixture --cv {params.out_prefix}.bed $i 1>{params.out_prefix}.$i.log\" 2>&1; done | parallel --progress -j {threads} ); "
		"min_CV=99999999999999; "
		"for i in {{{params.Kmin}..{params.Kmax}}}; do"
		"  grep \"CV\" {params.out_dir}/{params.out_prefix}.$i.log | awk '{{print $4}}' > {params.out_dir}/{params.out_prefix}.$i.CV;"
		"  CV=$(cat {params.out_dir}/{params.out_prefix}.$i.CV);"
		# Multiple by 1 million to get from float to int (unix if cant do floating point math)
		"  CV=$(printf %.0f $(echo \"$CV * 1000000\" | bc));"
		# Check if smaller then min and not zero (zero could be due to "-nan" in CV file - possible error/bug in admixture)
		"  if [[ $CV -lt $min_CV && $CV -ne 0 ]]; then"
		"    min_CV=$CV;"
		"    cp {params.out_dir}/{params.out_prefix}.$i.Q {output.best};"
		"  fi; "
		"done"
		") 1>{log} 2>&1"


rule format_results_plink_PCA:
	input:
		eigenval=rules.relatedness_plink_PCA.output.eigenval,
		eigenvec=rules.relatedness_plink_PCA.output.eigenvec,
	output:
		eigenval="results/{project}/final/PCA.eigenval",
		eigenvec="results/{project}/final/PCA.eigenvec",
	log:
		"results/logs/{project}/relatedness/format_plink_PCA_results.log",
	conda:
		"../envs/bash.yaml"
	shell:
		"("
		"cp {input.eigenval} {output.eigenval}; "
		"cp {input.eigenvec} {output.eigenvec}"
		") 1>{log} 2>&1"


rule format_results_plink_Admixture:
	input:
		rules.relatedness_plink_PCA.output.nosex,
		rules.relatedness_plink_Admixture.output.best,
	output:
		"results/{project}/final/Admixture.tsv",
	log:
		"results/logs/{project}/relatedness/format_plink_Admixture_results.log",
	conda:
		"../envs/bash.yaml"
	script:
		"../scripts/format_plink_Admixture_results.sh"


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


