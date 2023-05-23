#!/usr/bin/env python
import sys
import os
import subprocess
import click
import psutil, math
from snakemake.io import load_configfile
from snakemake.utils import validate


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

## Return max system memory or check if user specified meory is < system max
def get_system_max_mem():
	system_max_mem = math.floor(psutil.virtual_memory().total / (1024**3))
	system_max_mem = int(system_max_mem)
	return(system_max_mem)

def check_user_max_mem(user_max_mem):
	"""Check that user specified mem is <= system max mem.
	If user mem > system max mem return error.
	"""
	system_max_mem = get_system_max_mem()
	
	if (user_max_mem > system_max_mem):
		print("ERROR: You specified {user_max_mem}GB as maximum memory, but your system only has {system_max_mem}GB".format(
			user_max_mem=user_max_mem,
			system_max_mem=system_max_mem,
			)
		)
		sys.exit(1)



##
## Pass command line arguments.
##
@click.command(context_settings=dict(
	ignore_unknown_options=True, 
	show_default=True,
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
	"--max-memory",
	type=int,
	required=False,
	default=get_system_max_mem(),
	help="Max memory (RAM) allocated to the workflow. Snakemake will make sure jobs with memory allocations do not excede this limit.",
)
@click.option(
	"--max-downloads",
	required=False,
	default=6,
	help="Max number of parallel fasterq-dump jobs to run at the same time. Lets the user limit the amount that is being downloaded at the same time.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def workflow(module, configfile, cores, max_downloads, max_memory, snakemake_args):
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
	
	### Check max mem is <= system mem
	check_user_max_mem(max_memory)
	
	### Run Snakemake command
	run_cmd((
		"snakemake"
		" --config 'module={module}'"
		" --configfile '{configfile}'"
		" --cores {cores}"
		" --resources"
		" max_downloads={max_downloads}"
		" mem_gb={max_memory}"
		" {snakemake_args}"
	).format(
		module=module,
		configfile=configfile,
		cores=cores,
		max_downloads=max_downloads,
		max_memory=max_memory,
		snakemake_args=" ".join(SNAKEMAKE_REQUIRED_PARAMS + list(snakemake_args)),
	))
	
	#### Build report '.zip' file if we have finished running the workflow AND we havent previously made the report
	conf = load_configfile(configfile)
	project_name = conf["project_name"]
	
	done_workflow_path = "results/{project_name}/{module_name}.done".format(
		project_name=project_name, 
		module_name=module,
	)
	report_path   = "results/{project_name}/report.zip".format(
		project_name=project_name
	)
	
	# Check if report doesnt exist or if the workflow was run more reciently then the report was created.
	done_workflow_time = 0
	report_time = 0
	if os.path.isfile(done_workflow_path):
		done_workflow_time = os.path.getmtime(done_workflow_path)
	if os.path.isfile(report_path):
		report_time = os.path.getmtime(report_path)
	
	if done_workflow_time > report_time:
		run_cmd((
			"snakemake"
			" --config 'module={module}'"
			" --configfile '{configfile}'"
			" --report {report_path}"
			" {snakemake_args}"
		).format(
			module=module,
			configfile=configfile,
			report_path=report_path,
			snakemake_args=" ".join(SNAKEMAKE_REQUIRED_PARAMS + list(snakemake_args)),
		))



if __name__ == "__main__":
	workflow()
