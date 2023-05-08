# Genotype Samples
Snakemake workflow to genotype any type of NGS samples.


## Setup and installation
To install the workflow simply pull it from githib.
```bash
git clone https://github.com/TimothyStephens/Genotype_Samples.git
```

Then setup the conda environment which we will use to run the snakemake workflow.
```bash
conda env create -f environment.yaml
conda activate Genotype_Samples
```


## Workflow config
To run the workflow you will need to have two config/metadata files prepared.

The first config/metadata file is the `config/config.*.yaml` file.
 - The `config/config.Genotyping.yaml` file is a template for the genotyping workflow
 - The `config/config.crossMapping.yaml` file is a template for the cross_mapping workflow
When you want to run one of thse workflows make a copy of the relevant config file and fill it out.
The main params that need to be changed are:
 - `project_name`  The name of this project. A directory with the same name will be created in `results/` and is where the final analysis will be stored.
 - `samples`:      The path to the `samples.tsv` file to use for analysis.
 - `ref_genomes`:  The name and location of the reference genome to use for analysis. 

The second config/metadata file is the `config/samples.tsv` file.
This file lists the names and locations of all samples that will be considered by the workflow. 
 - `sample_id`  The ID of the sample
 - `unit`       A number used to differentiate samples with multiple input read sets. If you only have one set of reads per sample then keep this number as `1`. If you have multiple seqs of reads per sample then give each set a different number. Reads from the same sample_id but different "units" will be QC'ed, trimmed, and aligned separatly, but combine after that for final variant calling, genotyping, and ploidy analysis. 
 - `lib_type`   Lets you specify if the reads are from "dna" or "rna" sequencing. 
 - `fq1`        Location (absolute path) of the 1st read file OR the SRA ID of the sample which you wish the workflow to download and process for you.
 - `fq2`        Location (absolute path) of the 2nd read file (if present) OR the SRA ID of the sample which you wish the workflow to download and process for you. If you have samples that are single end, omit this column for those samples.

NOTE:
 - For SRA ID samples. You need to tell the workflow to expect the downloaded files to be single or paired-end. To specify paired-end you simply need to list the SRA ID in the `fq1` and `fq2` columns, to specify single-end specify the SRA ID in just the `fq1` column and omit the `fq2` column.
 - You can run each workflow multiple times with different subsets of samples. Just change the `project_name` in the config file and snakemake will do the rest. 
 - If you wish to add more samples to an existing project, simply add them to the config file and let snakemake do the rest.
 - Default config values are stored in `workflow/config.default.yaml`. If you want to modify them globally then do so in the `config.default.yaml` file, if you want to do it for a single project, simply add them to your project config file. The project config file provided by the user will overwrite any params in the global default file, so it is OK for a param to be set in both.
 - In the `samples.tsv` file, the `unit` column allows samples split across multiple `fastq` files to be combined.


## Running the workflow
This workflow is broken down into two parts that can be run separatly but using the same underlying datasets.
The `Genotype_Samples.py` script handels the high level running and organization of the workflow.

### cross_mapping
Compare samples against multiple genomes and report mapping stats (% mapped reads) as a way of checking/determining the "species" that a sample is from/has highest similarity to.
```bash
conda activate Genotype_Samples
./Genotype_Samples.py --module cross_mapping --configfile config/config.crossMapping.yaml
```
Once the workflow has finished you can run the same command with the `--report report.zip` flag to produce a nice HTML report.
```bash
./Genotype_Samples.py --module cross_mapping --configfile config/config.crossMapping.yaml --report report.zip
```
The HTML report will be contained within the `report.zip` file (which file can be named whatever you want in the above command so you dont overwrite reports from other projects).


### genotyping
Full genotyping workflow. Call variants, calculate relatedness stats, estimate ploidy, and generate QC and results plots.
```bash
conda activate Genotype_Samples
./Genotype_Samples.py --module genotyping --configfile config/config.Genotyping.yaml
```
To generate the HTML report run:
```bash
./Genotype_Samples.py --module genotyping --configfile config/config.Genotyping.yaml --report report.zip
```



## Run tests
Run test analyses using the following commands.
```bash
conda activate Genotype_Samples

# Cross-mapping
./Genotype_Samples.py --module cross_mapping --configfile tests/config/config.crossMapping_small.yaml
./Genotype_Samples.py --module cross_mapping --configfile tests/config/config.crossMapping_big.yaml

# Genotyping
./Genotype_Samples.py --module genotyping --configfile tests/config/config.Genotyping_small.yaml
./Genotype_Samples.py --module genotyping --configfile tests/config/config.Genotyping_big.yaml

# K-mer analysis
./Genotype_Samples.py --module kmer_analysis --configfile tests/config/config.KmerAnalysis_small.yaml
```
To generate HTML reports for each of the above commands, add `--report report.zip` to the end of each command.



## Notes
If you want to stop the scheduling of new jobs and wait for all running jobs to be finished, you can send a TERM signal, e.g., via
```bash
killall -TERM snakemake
```



## Results
Results files produced by this workflow will be in `results/project_name/final/`, where `project_name` is whatever you specified in the config file.

### `Annotations.tsv`
The ploidy and "group" (determined by vcf_clone_detect.py using the user specified threshold) of the input samples.

### `NgsRelate.tsv`
Relatedness results from the `NgsRelate` program. 

### `nQuire.tsv`
Ploidy estimated by `nQuire`. Estimates are for both raw (normal) and denoised biallelic sites.

### `PCAngsd.Admixture.tsv`
Admixture results produced by `PCAngsd` - number of ancestral populations chosen by `PCAngsd`.

### `PCAngsd.IndAlleleFreq.tsv`
PCA results produced by  `PCAngsd` with individual allele frequencies estimated.

### `PCAngsd.WithOutIndAlleleFreq.tsv`
PCA results produced by  `PCAngsd` without individual allele frequencies estimated.

### `vcf_clone_detect.tsv`
Matrix of the percent shared variants between each pairwise combination of samples, calculated by `vcf_clone_detect.py`.

### `vcftools_relatedness2.tsv`
Matrix of the relatedness values between each pairwise combination of samples, calculated by `vcftools`.


# Credits
This workflow is based on the [deepvariant workflow](https://github.com/nikostr/dna-seq-deepvariant-glnexus-variant-calling).
It also takes inspiration from the design of [Metagenome-Atlas](https://github.com/metagenome-atlas/atlas).

