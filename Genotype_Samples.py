#!/usr/bin/env python
import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import subprocess

__version__ = "0.0.1"



##
## Pass command line arguments.
##
DESCRIPTION = '''

Genome type samples

A snakemake workflow for exploring genotype + ploidy of DNA or RNA samples using a reference genome.

'''
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
subparsers = parser.add_subparsers(dest='command', required=True)



##
## Parser for the full typing (genotyping + ploidy) workflow
##
FULL_TYPING_WORKFLOW_DESCRIPTION = '''

Run full typing (genotyping + ploidy) workflow.

To pass options to snakemake put them in quotes at the end of this command.

'''
parser_full_typing_workflow = subparsers.add_parser('full_workflow', 
	help='Run full typing (genotyping + ploidy) workflow',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=FULL_TYPING_WORKFLOW_DESCRIPTION,
)
parser_full_typing_workflow.add_argument('--configfile', metavar='config/config.yaml', 
	required=True, type=str,
	help='Config file to use for analysis'
)
parser_full_typing_workflow.add_argument('-c', '--cores', default='all',
	required=False, type=str,
	help='Snakemake: Use at most N CPU cores/jobs in parallel (default: %(default)s)'
)
parser_full_typing_workflow.add_argument('snakemake_args', 
	nargs=argparse.REMAINDER
)



##
## Parser for the cross mapper workflow
##
CROSS_MAPPING_WORKFLOW_DESCRIPTION = '''

Run cross mapping workflow.

To pass options to snakemake put them in quotes at the end of this command.

'''
parser_cross_mapping_workflow = subparsers.add_parser('cross_mapping',
	help='Run cross mapping workflow',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=CROSS_MAPPING_WORKFLOW_DESCRIPTION,
)
parser_cross_mapping_workflow.add_argument('--configfile', metavar='config/config.yaml',
	required=True, type=str,
	help='Config file to use for analysis'
)
parser_cross_mapping_workflow.add_argument('-c', '--cores', default='all',
	required=False, type=str,
	help='Snakemake: Use at most N CPU cores/jobs in parallel (default: %(default)s)'
)
parser_cross_mapping_workflow.add_argument('snakemake_args',
	nargs=argparse.REMAINDER
)



##
## Parser for snakemake report
##
REPORT_DESCRIPTION = '''

Build report zip file.

To pass options to snakemake put them in quotes at the end of this command.

'''
parser_report_workflow = subparsers.add_parser('report',
	help='Build snakemake report',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=REPORT_DESCRIPTION,
)
parser_report_workflow.add_argument('--configfile', metavar='config/config.yaml',
	required=True, type=str,
	help='Config file to use for analysis'
)
parser_report_workflow.add_argument('--report',
	required=False, type=str, default="report.zip",
	help='Output report file name (default: %(default)s)'
)
parser_report_workflow.add_argument('snakemake_args',
        nargs=argparse.REMAINDER
)


##
## Parse all arguments.
##
args = parser.parse_args()

if args.command == "report":
	cmd = (
	"snakemake"
	" --configfile '{configfile}'"
	" --report {report}"
	" {snakemake_args}"
	).format(
		configfile=args.configfile,
		report=args.report,
		snakemake_args=" ".join(args.snakemake_args),
	)
else:
	cmd = (
		"snakemake --use-conda --use-singularity --keep-going"
		" --configfile '{configfile}'"
		" --cores {cores}"
		" {target_rule}"
		" {snakemake_args}"
	).format(
		configfile=args.configfile,
		cores=args.cores,
		target_rule=args.command,
		snakemake_args=" ".join(args.snakemake_args),
	)

print(cmd)
try:
	subprocess.check_call(cmd, shell=True)
except subprocess.CalledProcessError as e:
	# removes the traceback
	print("ERROR encountered while running snakemake:")
	print(e)
	exit(1)


