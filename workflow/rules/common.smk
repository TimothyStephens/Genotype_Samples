from snakemake.utils import validate
import pandas as pd
import numpy as np
import os


#######################
##### Load config #####
#######################
validate(config, schema="../schemas/config.schema.yaml")
PROJECT=config['project_name']


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


