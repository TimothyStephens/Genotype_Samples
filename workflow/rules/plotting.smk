

rule plotting_ANGSD_results:
	input:
		samples=rules.combine_genotyping_results.output.samples,
		color_list=rules.format_annotations.output.color_list,
		admixture=rules.format_results_PCAngsd_Admixture.output,
		allelFreqs=rules.format_results_PCAngsd_IndAlleleFreq.output,
	output:
		html=report(
			"results/{project}/final/ANGSD_results.html",
			#caption="../report/multiqc_calls.rst",
			subcategory="Result plots",
			labels={"Results": "ANGSD results"},
		),
		rmd="results/{project}/final/ANGSD_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_ANGSD_results.log",
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
		samples=rules.combine_genotyping_results.output.samples,
		color_list=rules.format_annotations.output.color_list,
		matrix=rules.format_results_vcf_clone_detect_counts.output,
	output:
		html=report(
			"results/{project}/final/vcf_clone_detect_results.html",
			#caption="../report/multiqc_calls.rst",
			subcategory="Result plots",
			labels={"Results": "vcf_clone_detect.py results"},
		),
		rmd="results/{project}/final/vcf_clone_detect_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_vcf_clone_detect_results.log",
	params:
		threshold=config["relatedness_vcf_clone_detect"]["threshold"],
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; sed -e 's/<<<similarity.threshold>>>/{params.threshold}/' workflow/scripts/plot_vcf_clone_detect_results.Rmd > {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"



rule plotting_vcftools_relatedness2_results:
	input:
		samples=rules.combine_genotyping_results.output.samples,
		color_list=rules.format_annotations.output.color_list,
		matrix=rules.format_results_vcftools_relatedness2.output,
	output:
		html=report(
			"results/{project}/final/vcftools_relatedness2_results.html",
			#caption="../report/multiqc_calls.rst",
			subcategory="Result plots",
			labels={"Results": "vcftools --relatedness2 results"},
		),
		rmd="results/{project}/final/vcftools_relatedness2_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_vcftools_relatedness2_results.log",
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


rule plotting_nQuire_delta_loglikelihood_results:
	input:
		rules.format_nQuire_results.output,
	output:
		html=report(
			"results/{project}/final/nQuire_delta_loglikelihood_results.html",
			#caption="../report/multiqc_calls.rst",
			subcategory="Result plots",
			labels={"Results": "nQuire delta loglikelihood results"},
		),
		rmd="results/{project}/final/nQuire_delta_loglikelihood_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_nQuire_delta_loglikelihood_results.log",
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; cp workflow/scripts/plot_nQuire_delta_loglikelihood_results.Rmd {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"


rule plotting_nQuire_coverage_results:
	input:
		rules.ploidy_nQuire_overage_file_list.output,
	output:
		html=report(
			"results/{project}/final/nQuire_sites_coverage_results_xLim_0.1-0.9.pdf",
			#caption="../report/multiqc_calls.rst",
			subcategory="Result plots",
			labels={"Results": "nQuire sites coverage"},
		),
		rmd="results/{project}/final/nQuire_sites_coverage_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_nQuire_sites_coverage_results.log",
	params:
		ref_name=list(config["ref_genomes"].keys())[0],
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; sed "
		"    -e 's@<<<files2plot.prefix>>>@../../ploidy/{params.ref_name}/@'"
		"   workflow/scripts/plot_nQuire_sites_coverage_results.Rmd > {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"


rule plotting_cross_mapping_results:
	input:
		matrix=rules.format_crossMapping_results.output,
		samples=rules.format_annotations.output.samples,
		color_list=rules.format_annotations.output.color_list,
	output:
		html=report(
			"results/{project}/final/cross_mapping_results.html",
			#caption="../report/multiqc_calls.rst",
			subcategory="Result plots",
			labels={"Results": "Cross mapping results"},
		),
		rmd="results/{project}/final/cross_mapping_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_cross_mapping_results.log",
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; cp workflow/scripts/plot_cross_mapping_results.Rmd {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"


rule plotting_kmer_analysis_results:
	input:
		genomescape2_plots_linear=rules.kmer_analysis_genomescope2.output.plot_linear,
		genomescape2_plots_log=rules.kmer_analysis_genomescope2.output.plot_log,
		genomescape2_plots_transformed_linear=rules.kmer_analysis_genomescope2.output.plot_transformed_linear,
		genomescape2_plots_transformed_log=rules.kmer_analysis_genomescope2.output.plot_transformed_log,
		smudgeplot_plot=rules.kmer_analysis_smudgeplot.output.plot,
		smudgeplot_plot_log=rules.kmer_analysis_smudgeplot.output.plot_log,
	output:
		html=report(
			"results/{project}/kmer_analysis/{sample}.kmer_analysis_results.html",
			#caption="../report/multiqc_calls.rst",
			subcategory="K-mer analysis results",
			labels={"Sample ID": "{sample}"},
		),
		rmd="results/{project}/kmer_analysis/{sample}.kmer_analysis_results.Rmd",
	log:
		"results/logs/{project}/plotting/plot_kmer_analysis_results.{sample}.log",
	conda:
		"../envs/R.yaml"
	shell:
		"("
		"  export PATH=\"$CONDA_PREFIX/bin:$PATH\""
		"; export R_LIB=\"$CONDA_PREFIX/lib/R/library\""
		"; sed "
		"    -e 's@<<<genomescape2_plots_linear>>>@{input.genomescape2_plots_linear}@'"
		"    -e 's@<<<genomescape2_plots_log>>>@{input.genomescape2_plots_log}@'"
		"    -e 's@<<<genomescape2_plots_transformed_linear>>>@{input.genomescape2_plots_transformed_linear}@'"
		"    -e 's@<<<genomescape2_plots_transformed_log>>>@{input.genomescape2_plots_transformed_log}@'"
		"    -e 's@<<<smudgeplot_plot>>>@{input.smudgeplot_plot}@'"
		"    -e 's@<<<smudgeplot_plot_log>>>@{input.smudgeplot_plot_log}@'"
		"   workflow/scripts/plot_kmer_analysis.Rmd > {output.rmd}"
		"; Rscript -e \"rmarkdown::render('{output.rmd}')\""
		")"
		" 1>{log} 2>&1"


