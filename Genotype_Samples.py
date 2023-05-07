#!/usr/bin/env python
import sys
import os
import subprocess
import click


__version__ = "0.0.1"


## Params to be used by all snakemake commands.
SNAKEMAKE_REQUIRED_PARAMS = ["--use-conda", "--use-singularity", "--keep-going", "--printshellcmds"]

## Helper script to run snakemake commands.
def run_cmd(cmd):
	print("Running snakemake command:\n", cmd)
	try:
		subprocess.check_call(cmd, shell=True)
	except subprocess.CalledProcessError as e:
		# removes the traceback
		print("ERROR encountered while running snakemake:")
		print(e)
		exit(1)


##
## Pass command line arguments.
##
@click.command(context_settings=dict(
	ignore_unknown_options=True, 
	help_option_names=["-h", "--help"],
))
@click.version_option(__version__)
@click.option(
	"--module",
	required=True,
	type=click.Choice(['genotyping', 'cross_mapping', 'kmer_analysis']),
	help="Module/workflow to run",
)
@click.option(
	"--configfile",
	required=True,
	type=click.Path(exists=True, resolve_path=True),
	help="configfile",
)
@click.option(
	"-c",
	"--cores",
	required=False,
	default="all",
	help="Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the number of available CPU cores.",
)
@click.option(
	"--max-downloads",
	required=False,
	default=6,
	help="Max number of parallel fasterq-dump jobs to run at the same time. Lets the user limit the amount that is being downloaded at the same time.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def workflow(module, configfile, cores, max_downloads, snakemake_args):
	"""
	\b
	##########################
	#### Genotype samples ####
	##########################
	
	A snakemake workflow for exploring genotype + ploidy of DNA or RNA samples using a reference genome.
	
	\b
	Modules:
	  genotyping		Exploring genotype + ploidy of samples
	  cross_mapping		Cross-mapping of samples to multiple genomes
	  kmer_analysis		K-mer ploidy analysis
	
	##########################
	"""
	run_cmd((
		"snakemake"
		" --config 'module={module}'"
		" --configfile '{configfile}'"
		" --cores {cores}"
		" --resources max_downloads={max_downloads}"
		" {snakemake_args}"
	).format(
		module=module,
		configfile=configfile,
		cores=cores,
		max_downloads=max_downloads,
		snakemake_args=" ".join(SNAKEMAKE_REQUIRED_PARAMS + list(snakemake_args)),
	))



if __name__ == "__main__":
	workflow()
