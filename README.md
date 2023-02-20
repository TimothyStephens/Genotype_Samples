# Genotype_Samples
Snakemake workflow to genotype any type of NGS sample

Based on https://github.com/nikostr/dna-seq-deepvariant-glnexus-variant-calling with significant modifications.

mamba create -n snakemake snakemake mamba singularity
conda activate snakemake

snakemake --cores all --use-conda --use-singularity --keep-going 
snakemake --report report.zip


