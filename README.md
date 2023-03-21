# Genotype Samples
Snakemake workflow to genotype any type of NGS samples.

Based on https://github.com/nikostr/dna-seq-deepvariant-glnexus-variant-calling with significant modifications.



## Setup environment

Setup conda environment which we will use to run the snakemake workflow.

```bash
mamba create -n snakemake snakemake=7.22.0 singularity=3.8.6 mamba
conda activate snakemake
```



## Workflow config

Setup config files.

### `config/config.yaml`

Workflow config file. The only parameter which needs to be changed is the reference genome that is being used.

### `config/samples.tsv`

File listing the samples to analyze. 

NOTE: the `unit` column allows samples split across multiple `fastq` files to be combined. 

### `config/joint_calling_groups.tsv`

File listing samples to jointly genotype. Leave blank (i.e., just column headers) to genotype individually. 



## Run workflow

Use the following commands to run the workflow. 

```bash
snakemake --cores all --use-conda --use-singularity --keep-going 
```

Create run report.

```bash
snakemake --report report.zip
```



## Results

Results files produced by this workflow.

### `results/final/NgsRelate.tsv`

Relatedness results from the `NgsRelate` program. 

### `results/final/nQuire.tsv`

Ploidy estimated by `nQuire`. Estimates are for both raw (normal) and denoised biallelic sites.

### `results/final/PCAngsd.Admixture.tsv`

Admixture results produced by `PCAngsd` - number of ancestral populations chosen by `PCAngsd`.

### `results/final/PCAngsd.IndAlleleFreq.tsv`

PCA results produced by  `PCAngsd` with individual allele frequencies estimated.

### `results/final/PCAngsd.WithOutIndAlleleFreq.tsv`

PCA results produced by  `PCAngsd` without individual allele frequencies estimated.

### `results/final/vcf_clone_detect.tsv`

Matrix of the percent shared variants between each pairwise combination of samples, calculated by `vcf_clone_detect.py`.

### `results/final/vcftools_relatedness2.tsv`

Matrix of the relatedness values between each pairwise combination of samples, calculated by `vcftools`.



`Rmarkdown` scripts for visualization of the `PCAngsd`, `vcftools_relatedness2`, and `vcf_clone_detect` results are in `workflow/scripts/R/`.


