from snakemake.utils import validate
import pandas as pd
import numpy as np
import os


#######################
##### Load config #####
#######################
validate(config, schema="../schemas/config.schema.yaml")


#############################
##### Load sample sheet #####
#############################
samples = pd.read_table(config["samples"], dtype=str, comment='#')
validate(samples, schema="../schemas/samples.schema.yaml")

# Set index of dataframe
samples = samples.set_index(
	["sample_id", "unit", "lib_type"], drop=False
)
# enforce str in index
samples.index = samples.index.set_levels(
	[i.astype(str) for i in samples.index.levels]
)


############################
##### Helper Functions #####
############################

def expand_raw_fastqc_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2):
			out.append("{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit))
			out.append("{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit))
		if pd.isnull(row.fq2):
			out.append("{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit))
	return out


def expand_fastq_paths():
	out = []
	for i, row in samples.iterrows():
		if row.lib_type == "dna-pe" or row.lib_type == "rna-pe": ## DNA or RNA pe
			out.append("{lib_type}/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("{lib_type}/{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		else: ## DNA and RNA se + DNA and RNA long
			out.append("{lib_type}/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


def expand_sample_paths():
	out = []
	for i, row in samples.iterrows():
		out.append("{lib_type}/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


