

rule plotting_ANGSD_results:
	input:
		annot=rules.format_annotations.output,
		admixture=rules.format_results_PCAngsd_Admixture.output,
		allelFreqs=rules.format_results_PCAngsd_IndAlleleFreq.output,
	output:
		html=report(
			"results/"+PROJECT+"/final/ANGSD_results.html",
			#caption="../report/multiqc_calls.rst",
			category=PROJECT,
			subcategory="Result Plots",
			labels={"Results": "ANGSD results"},
		),
		rmd="results/"+PROJECT+"/final/ANGSD_results.Rmd",
	log:
		"results/"+PROJECT+"/log/relatedness/plot_ANGSD_results.log",
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; cp workflow/scripts/plot_ANGSD_results.Rmd {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"


rule plotting_vcf_clone_detect_results:
	input:
		annot=rules.format_annotations.output,
		matrix=rules.format_results_vcf_clone_detect_counts.output,
	output:
		html=report(
			"results/"+PROJECT+"/final/vcf_clone_detect_results.html",
			#caption="../report/multiqc_calls.rst",
			category=PROJECT,
			subcategory="Result plots",
			labels={"Results": "'vcf_clone_detect.py' results"},
		),
		rmd="results/"+PROJECT+"/final/vcf_clone_detect_results.Rmd",
	log:
		"results/"+PROJECT+"/log/relatedness/plot_vcf_clone_detect_results.log",
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; cp workflow/scripts/plot_vcf_clone_detect_results.Rmd {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"



rule plotting_vcftools_relatedness2_results:
	input:
		annot=rules.format_annotations.output,
		matrix=rules.relatedness_vcftools_relatedness2.output,
	output:
		html=report(
			"results/"+PROJECT+"/final/vcftools_relatedness2_results.html",
			#caption="../report/multiqc_calls.rst",
			category=PROJECT,
			subcategory="Result plots",
			labels={"Results": "'vcftools --relatedness2' results"},
		),
		rmd="results/"+PROJECT+"/final/vcftools_relatedness2_results.Rmd",
	log:
		"results/"+PROJECT+"/log/relatedness/plot_vcftools_relatedness2_results.log",
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; cp workflow/scripts/plot_vcftools_relatedness2_results.Rmd {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"


