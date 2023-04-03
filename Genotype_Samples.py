#!/usr/bin/env python
import sys
import os
from argparse import RawTextHelpFormatter
import subprocess
import click

__version__ = "0.0.1"


def run_cmd(cmd):
	print(cmd)
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
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
	"""
	Genotype samples
	A snakemake workflow for exploring genotype + ploidy of DNA or RNA samples using a reference genome.
	"""

@cli.command(
	"genotyping",
	context_settings=dict(ignore_unknown_options=True),
	short_help="Run full typing (genotyping + ploidy) workflow",
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
def run_genotyping(configfile, cores, max_downloads, snakemake_args):
	run_cmd((
		"snakemake --use-conda --use-singularity --keep-going"
		" --snakefile '{snakefile}'"
		" --configfile '{configfile}'"
		" --cores {cores}"
		" --resources max_downloads={max_downloads}"
		" {snakemake_args}"
	).format(
		snakefile="workflow/Snakefile_genotyping",
		configfile=configfile,
		cores=cores,
		max_downloads=max_downloads,
		snakemake_args=" ".join(snakemake_args),
	))


@cli.command(
	"cross_mapping",
	context_settings=dict(ignore_unknown_options=True),
	short_help="Run cross mapping workflow",
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
def run_cross_mapping(configfile, cores, max_downloads, snakemake_args):
	run_cmd((
		"snakemake --use-conda --use-singularity --keep-going"
		" --snakefile '{snakefile}'"
		" --configfile '{configfile}'"
		" --cores {cores}"
		" --resources max_downloads={max_downloads}"
		" {snakemake_args}"
	).format(
		snakefile="workflow/Snakefile_cross_mapping",
		configfile=configfile,
		cores=cores,
		max_downloads=max_downloads,
		snakemake_args=" ".join(snakemake_args),
	))


if __name__ == "__main__":
	cli()
