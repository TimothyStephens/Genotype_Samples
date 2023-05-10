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
samples = pd.read_table(config["samples"], dtype=str).set_index(
	["sample_id", "unit", "lib_type"], drop=False
)
samples.index = samples.index.set_levels(
	[i.astype(str) for i in samples.index.levels]
)  # enforce str in index
validate(samples, schema="../schemas/samples.schema.yaml")


############################
##### Helper Functions #####
############################

def expand_raw_fastqc_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2):
			out.append("{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		if pd.isnull(row.fq2):
			out.append("{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


def expand_fastq_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2): ## DNA pe
			out.append("{lib_type}/pe/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
			out.append("{lib_type}/pe/{sample}-{unit}.2".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		else:
			if row.lib_type == "dna-long": ## DNA long
				out.append("{lib_type}/long/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type="dna"))
			elif row.lib_type == "rna-long": ## RNA long
				out.append("{lib_type}/long/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type="rna"))
			else: ## DNA se
				out.append("{lib_type}/se/{sample}-{unit}.1".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


def expand_sample_paths():
	out = []
	for i, row in samples.iterrows():
		if pd.notnull(row.fq2): ## DNA pe
			out.append("{lib_type}/pe/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
		else:
			if row.lib_type == "dna-long": ## DNA long
				out.append("{lib_type}/long/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type="dna"))
			elif row.lib_type == "rna-long": ## RNA long
				out.append("{lib_type}/long/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type="rna"))
			else: ## DNA se
				out.append("{lib_type}/se/{sample}-{unit}".format(sample=row.sample_id, unit=row.unit, lib_type=row.lib_type))
	return out


