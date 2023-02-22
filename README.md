# Genotype_Samples
Snakemake workflow to genotype any type of NGS sample

Based on https://github.com/nikostr/dna-seq-deepvariant-glnexus-variant-calling with significant modifications.

mamba create -n snakemake snakemake=7.22.0 singularity=3.8.6 mamba
conda activate snakemake

snakemake --cores all --use-conda --use-singularity --keep-going 
snakemake --report report.zip


