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
 - Changed the "group" identification step to now use a network connected 
   component identification algorithm. Is more robust when we have two sets
   of weakly related samples that have similarity above the threshold but the 
   previous algorithm would split into two groups.

Attempts to identify groups of clones in a dataset. The script
(1) LOAD pairwise comparisons (allelic similarity) produced by `vcf_clone_detect-count.py` 
	for all individuals,
(2) produces a histogram of genetic similarities,
(3) lists the highest matches to assess for a potential clonal threshold,
(4) clusters the groups of clones based on a particular threshold (supplied or
	roughly inferred; group is defined as any sample with similairty to another
        sample in the group > user threshold), and
(5) lists the clonal individuals that can be removed from the dataset
	(so that one individual with the least amount of missing data remains).

e.g. `python3 vcf_clone_detect-group.py --input compare_file.csv --output sample_groups.tsv --threshold 95`

"""
import sys
import argparse
import operator
import itertools
import math
import numpy as np
from scipy import sparse
from sknetwork.data import from_edge_list
from sknetwork.topology import get_connected_components

__author__ = 'Pim Bongaerts'
__copyright__ = 'Copyright (C) 2016 Pim Bongaerts'
__license__ = 'GPL'


IND_STR_LEN = 50
MATCHES_ADDITIONAL_ROWS = 5
HIST_RANGE = range(1, 101)
DEF_THRESHOLD = 95.0

C_IND1 = 'ind1'
C_IND2 = 'ind2'
C_IND1_SNPS = 'ind1_snps'
C_IND2_SNPS = 'ind2_snps'
C_BOTH_SNPS = 'both_snps'
C_MATCH = 'match'
C_MATCH_PERC = 'match_perc'

COMPARISONS_DTYPES = [(C_IND1, np.str_, IND_STR_LEN),
	(C_IND2, np.str_, IND_STR_LEN),
	(C_IND1_SNPS, int),
	(C_IND2_SNPS, int),
	(C_BOTH_SNPS, int),
	(C_MATCH, float),
	(C_MATCH_PERC, float)]


class CloneGroup(object):

	def __init__(self, row):
		self.indivs = set([str(row[C_IND1]), str(row[C_IND2])])
		self.min_sim_score = self.max_sim_score = float(row[C_MATCH_PERC])
		if int(row[C_IND1_SNPS]) >= int(row[C_IND2_SNPS]):
			self.best_indiv = str(row[C_IND1])
			self.best_indiv_snps = int(row[C_IND1_SNPS])
		else:
			self.best_indiv = str(row[C_IND2])
			self.best_indiv_snps = int(row[C_IND2_SNPS])

	def add_clone_from_row(self, row):
		self.indivs.update([str(row[C_IND1]), str(row[C_IND2])])
		if float(row[C_MATCH_PERC]) < self.min_sim_score:
			self.min_sim_score = float(row[C_MATCH_PERC])
		if row[C_MATCH_PERC] > self.max_sim_score:
			self.max_sim_score = float(row[C_MATCH_PERC])
		if row[C_IND1_SNPS] >= self.best_indiv_snps:
			self.best_indiv = str(row[C_IND1])
			self.best_indiv_snps = int(row[C_IND1_SNPS])
		if row[C_IND2_SNPS] >= self.best_indiv_snps:
			self.best_indiv = str(row[C_IND2])
			self.best_indiv_snps = int(row[C_IND2_SNPS])

	def get_formatted_clone_info(self):
		if self.min_sim_score == self.max_sim_score:
			score_range = '{0}%'.format(self.min_sim_score)
		else:
			score_range = '{0}-{1}%'.format(self.min_sim_score, self.max_sim_score)
		info = '{0}\t({1})'.format(score_range, ', '.join(self.indivs))
		return info

	def get_samples_to_remove(self):
		return self.indivs - set([self.best_indiv])
	
	def get_samples(self):
		return self.indivs


def get_pairwise_comparisons_from_input_file(input_filename):
	""" Load pairwise comparisons between all individuals from input_file """
	comparisons = np.genfromtxt(input_filename, dtype=COMPARISONS_DTYPES, delimiter=',')
	comparisons[::-1].sort(order=C_MATCH_PERC)
	print('{0} comparisons loaded from `{1}`'.format(comparisons.size, input_filename))
	return comparisons


def output_ascii_hist(raw_values, bin_values):
	""" Plot a text-based histogram """
	values, bins = np.histogram(raw_values, bins=bin_values)
	output_lines = []
	graph_multiplier = 1
	lower_bound_flag = False
	lower_bound = 100 ##TGS
	upper_bound = 100 ##TGS
	previous_value = breakpoint = display_lines = 0
	for index, value in enumerate(values):
		if value > 0:
			if not lower_bound_flag:
				lower_bound_flag = True
				lower_bound = bins[index]
			else:
				upper_bound = bins[index]
		if lower_bound_flag:
			graph_bar = '*' * int(value * graph_multiplier)
			if len(graph_bar) > 70:
				graph_bar = graph_bar[:69] + '#'
			output_lines.append(
				'{:3d} {:7d} {:70s}'.format(bins[index], value, graph_bar[:70]))
	print('\n'.join(output_lines[:(upper_bound - lower_bound + 2)]))


def output_highest_matches(comparisons, threshold):
	""" Output list of highest matches """
	extra_rows = last_value = diff = highest_diff = 0
	highest_diff_max_perc = highest_diff_min_perc = 0
	output_lines = []
	if len(comparisons) > 0: ## TGS - Only run if not empty (prevents numpy error when iter empty arrays)
		for row in np.nditer(comparisons):
			# Only output matches above threshold (+ several additional rows)
			if threshold == 0 and row[C_MATCH_PERC] < DEF_THRESHOLD:
				break
			if threshold > 0 and row[C_MATCH_PERC] < threshold:
				extra_rows += 1
				#if extra_rows > MATCHES_ADDITIONAL_ROWS:
				#	break
			# Keep track of largest difference between sequential matches
			if last_value != 0:
				diff = round(last_value - row[C_MATCH_PERC], 2)
				if diff > highest_diff:
					highest_diff = diff
					highest_diff_min_perc = row[C_MATCH_PERC]
					highest_diff_max_perc = last_value
			
			output_lines.append((
				'{0}\t{1}\t{2} vs {3}\t{4}/{5}\t{6}\t{7}').format(
					row[C_MATCH_PERC],
					diff,
					row[C_IND1], 
					row[C_IND2],
					row[C_MATCH],
					row[C_BOTH_SNPS],
					row[C_IND1_SNPS],
					row[C_IND2_SNPS]
				)
			)
			last_value = row[C_MATCH_PERC]
	
	# Determine threshold value
	if threshold > 0:
		threshold_msg = 'Manual threshold'
	elif threshold == 0:
		threshold_msg = 'Potential threshold'
		if int(highest_diff_max_perc) > highest_diff_min_perc:
			threshold = float(int(highest_diff_max_perc))
		else:
			threshold = (highest_diff_max_perc - highest_diff_min_perc) / 2

	# Output list of matches with a break at the threshold
	for line in output_lines:
		if float(line.split()[0]) < threshold and threshold_msg:
			print('{0}\t{1} {2} {1}'.format(round(threshold, 2), '-' * 20, threshold_msg))
			threshold_msg = ''  # Use as flag so threshold only occurs once
		print(line)
	return threshold


def cluster_clones(comparisons, threshold):
	""" Cluster groups of clones together """
	clone_groups = []        # list of CloneGroup instances
	clone_indexes = {}       # clone_indexes[indiv] = (index for clone_groups)
	edge_list = []
	
	# End early if empty (prevents numpy error when iter empty arrays)
	if len(comparisons) == 0:
		return clone_groups
	
	for row in np.nditer(comparisons):
		# Only take pairs (edges) if they are above our chosen threshold.
		if row[C_MATCH_PERC] >= threshold:
			edge_list.append(row)
	
	# Edges to graph object
	graph = from_edge_list([ (str(row[C_IND1]), str(row[C_IND2])) for row in edge_list ])
	# Get connected components
	labels = get_connected_components(graph.adjacency)
	
	# Annotate each label with its component index
	for i, name in enumerate(graph.names):
		clone_indexes[name] = labels[i]
	
	# Setup empty list
	for i in range(max(labels)+1):
		clone_groups.append(None)
	
	# Add info for each group into CloneGroup instances
	for row in edge_list:
		ind1 = str(row[C_IND1])
		ind2 = str(row[C_IND2])
		ind1_idx = clone_indexes[ind1]
		ind2_idx = clone_indexes[ind2]
		if clone_groups[ind1_idx] is not None:
			clone_groups[ind1_idx].add_clone_from_row(row)
		else:
			clone_groups[ind1_idx] = CloneGroup(row)
	return clone_groups


def write_sample_groups(clone_groups, comparisons, output_filename):
	# Get list of all sample names.
	sample_ids = []
	for row in np.nditer(comparisons):
		sample_ids.append(str(row[C_IND1]))
		sample_ids.append(str(row[C_IND2]))
	sample_ids = list(set(sample_ids))
	
	# Write biggest to smallest group; also write ungrouped samples.
	sample_groups = sorted([c.get_samples() for c in clone_groups], key=lambda l: (len(l), l), reverse=True)
	clone_group = 1
	with open(output_filename, 'w') as fh:
		fh.write('sample_id\tgroup_id\n')
		for group in sample_groups:
			print(group)
			for sample in group:
				fh.write('{}\tGroup{}\n'.format(sample, clone_group))
				sample_ids.remove(sample)
			clone_group += 1
		for sample in sample_ids:
			fh.write('{}\tUngroup\n'.format(sample))
	

def main(input_filename, output_filename, threshold):
	print('### 1 - LOAD pairwise comparisons of all individuals')
	comparisons = get_pairwise_comparisons_from_input_file(input_filename)
	
	print('\n### 2 - Histogram (of pairwise genetic similarities)')
	output_ascii_hist(comparisons[C_MATCH_PERC], HIST_RANGE)
	
	print('\n### 3 - List of highest matches')
	threshold = output_highest_matches(comparisons, float(threshold))
	
	print('\n### 4 - Clonal groups (threshold: {0})'.format(threshold))
	clone_groups = cluster_clones(comparisons, float(threshold))
	for clone_group in clone_groups:
		print(clone_group.get_formatted_clone_info())
	
	print('\n### 5 - Individuals to remove from dataset (retaining indiv with least amount of missing data)')
	for clone_group in clone_groups:
		print('\n'.join(clone_group.get_samples_to_remove()))
	
	print('\n### 6 - Print groupings for all samples')
	write_sample_groups(clone_groups, comparisons, output_filename)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i', '--input', dest='input_filename',
		metavar='compare_file',
		required=True, type=str,
		help='input file (csv) with previously calculated pairwise comparisons (using the outputfile from `vcf_clone_detect-count.py`)')
	parser.add_argument('-o', '--output', dest='output_filename',
		metavar='compare_file',
		required=True, type=str,
		help='output file with groupings for all samples (tsv)' )
	parser.add_argument('-t', '--threshold', dest='threshold',
		metavar='threshold', default=95.0,
		required=False, type=float,
		help='manual similarity threshold (e.g. `95` means at least 95 percent allelic similarity for individuals to be considered clones)')
	args = parser.parse_args()
	main(args.input_filename, args.output_filename, args.threshold)
