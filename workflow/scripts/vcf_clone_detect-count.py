#!/usr/bin/env python
"""
##########################
#### vcf_clone_detect ####
##########################

Script (and all credit) originally by Pim Bongaerts: https://github.com/pimbongaerts/radseq

Modified to:
 - Fix minor bugs that arise then empty or sparse VCFs are given as input.
 - Allow for the SNP counting and clone detection to be performed as separate 
   scripts. This lets us rerun the clonal group detection using a different 
   cutoff as part of the snakemake worklfow. This functionality is already
   part of the script but having this as two separate scripts works better
   with the logic of snakemake.

Attempts to identify groups of clones in a dataset. The script
(1) conducts pairwise comparisons (allelic similarity) for all individuals in
    a `.vcf` file,

e.g. `python3 vcf_clone_detect-count.py.py --vcf SNPs.vcf --output compare_file.csv`

"""
import sys
import argparse
import operator
import itertools
import math
import numpy as np

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


HEADER_CHAR = "#"
HEADER_INDIVIDUALS = "#CHROM"
FIRST_GENOTYPE_COLUMN = 9
IND_STR_LEN = 50
POP_STR_LEN = 50
MATCHES_ADDITIONAL_ROWS = 5
HIST_RANGE = range(1, 101)
DEF_THRESHOLD = 85.0

C_IND1 = 'ind1'
C_IND2 = 'ind2'
C_IND1_SNPS = 'ind1_snps'
C_IND2_SNPS = 'ind2_snps'
C_BOTH_SNPS = 'both_snps'
C_MATCH = 'match'
C_MATCH_PERC = 'match_perc'
C_POP = 'pop'

COMPARISONS_DTYPES = [(C_IND1, np.str_, IND_STR_LEN),
	(C_IND2, np.str_, IND_STR_LEN),
	(C_IND1_SNPS, int),
	(C_IND2_SNPS, int),
	(C_BOTH_SNPS, int),
	(C_MATCH, float),
	(C_MATCH_PERC, float),
	(C_POP, np.str_, POP_STR_LEN)]

OUTPUT_FILE_DELIM = ','
OUTPUT_FILE_HEADER = OUTPUT_FILE_DELIM.join([C_IND1, C_IND2, C_IND1_SNPS,
	C_IND2_SNPS, C_BOTH_SNPS,
	C_MATCH, C_MATCH_PERC, C_POP])
OUTPUT_FILE_FORMAT = OUTPUT_FILE_DELIM.join(['%s', '%s', '%i', '%i', '%i', '%f', '%f', '%s'])


def get_snp_match(genotype1, genotype2):
	""" Get match value for two genotypes (one SNP) """
	if genotype1 == genotype2:
		match_score = 1
	elif genotype1[0] == genotype2[0] or genotype1[2] == genotype2[2] or \
			genotype1[0] == genotype2[2] or genotype1[2] == genotype2[0]:
		match_score = 0.5
	else:
		match_score = 0
	return match_score



def get_genotypes_from_vcf(vcf_filename):
	""" Read genotypes from vcf and store in dictionary """
	vcf_file = open(vcf_filename, 'r')
	individuals = {}
	genotypes = {}

	# Iterate through vcf file and store individual names and genotypes
	for line in vcf_file:
		cols = line.split()
		# Store name of individuals in temporary dictionary
		if line[:len(HEADER_INDIVIDUALS)] == HEADER_INDIVIDUALS:
			for x in range(FIRST_GENOTYPE_COLUMN, len(cols)):
				individuals[x] = cols[x]
		# Store genotype for each individual in dictionary
		elif line[0] != HEADER_CHAR:
			for x in range(FIRST_GENOTYPE_COLUMN, len(cols)):
				genotypes.setdefault(individuals[x], []).append(cols[x])
	vcf_file.close()
	return genotypes


def get_match_info_as_row(individual1, genotypes1,
						  individual2, genotypes2, indivs_pops):
	""" Get match value for two genotypes (one SNP) """
	ind1_snps = ind2_snps = both_snps = match = total_snps = 0
	for x in range(0, len(genotypes1)):
		genotype1 = genotypes1[x][:3]
		genotype2 = genotypes2[x][:3]
		if genotype1[0] != '.' and genotype2[0] != '.':  # both genotyped
			ind1_snps += 1
			ind2_snps += 1
			both_snps += 1
			match += get_snp_match(genotype1, genotype2)
		elif genotype2[0] == '.':  # genotype of ind2 missing
			ind1_snps += 1
		elif genotype1[0] == '.':  # genotype of ind1 missing
			ind2_snps += 1

	match_perc = round((match / both_snps) * 100, 2)

	if individual1 in indivs_pops and individual2 in indivs_pops:
		# Define whether comparison is within or between pops
		if indivs_pops[individual1] == indivs_pops[individual2]:
			pop_group = indivs_pops[individual1]
		else:
			pops = sorted([indivs_pops[individual1],
						   indivs_pops[individual2]])
			pop_group = '-'.join(pops)
	else:
		# Not applicable as indiv(s) not in popfile or popfile not provided
		pop_group = 'NA'

	return (str(individual1), str(individual2), int(ind1_snps),
			int(ind2_snps), int(both_snps), float(match),
			float(match_perc), str(pop_group))


def get_all_pairwise_comparisons(genotypes, indivs_pops):
	""" Conduct pairwise comparisons between all individuals """
	unique_pairs = int((math.pow(len(genotypes), 2) - len(genotypes)) / 2)
	comparisons = np.zeros(unique_pairs, dtype=COMPARISONS_DTYPES)
	
	index = 0
	for individual1, individual2 in itertools.combinations(genotypes, 2):
		comparisons[index] = get_match_info_as_row(individual1,
			genotypes[individual1],
			individual2,
			genotypes[individual2],
			indivs_pops)
		index += 1
	comparisons[::-1].sort(order=C_MATCH_PERC)
	print('{0} comparisons completed'.format(comparisons.size))
	return comparisons



def save_pairwise_comparisons_to_input_file(output_filename, comparisons):
	""" Load pairwise comparisons between all individuals to output_file """
	np.savetxt(output_filename, comparisons, fmt=OUTPUT_FILE_FORMAT,
		header=OUTPUT_FILE_HEADER)
	print('Comparisons outputted to file: `{0}`'.format(output_filename))




def main(vcf_filename, output_filename):
	print('### 1 - Pairwise comparisons of all individuals')
	indivs_pops = []
	
	# Genotypes
	genotypes = get_genotypes_from_vcf(vcf_filename)
	
	# Comparisons
	comparisons = get_all_pairwise_comparisons(genotypes, indivs_pops)
	
	# Save matrix data
	save_pairwise_comparisons_to_input_file(output_filename, comparisons)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,)
	parser.add_argument('-v', '--vcf', dest='vcf_filename', metavar='in.vcf',
		required=True, type=str,
		help='input file with SNP data (e.g., `*.vcf`)')
	parser.add_argument('-o', '--output', dest='output_filename',
		metavar='compare_file',
		required=True, type=str,
		help='output file (csv) for all pairwise comparisons')
	args = parser.parse_args()
	main(args.vcf_filename, args.output_filename)
