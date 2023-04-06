#!/usr/bin/env python2.7
#usage: 
# bwa mem -a assembly.fa reads_1.fastq reads_2.fastq | ./filter_reads_for_bin_reassembly.py original_bin_folder output_dir
import sys, os, gzip

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n'} 
def rev_comp(seq):
	rev_comp=""
	for n in seq:
		rev_comp+=complement[n]
	return rev_comp[::-1]

# load bin contigs
print "loading contig to bin mappings..."
contig_bins={}
contig_bins['*']="unalinged"
for bin_file in os.listdir(sys.argv[1]):
	if bin_file.endswith(".fa") or bin_file.endswith(".fasta"): 
		bin_name=".".join(bin_file.split("/")[-1].split(".")[:-1])
		for line in open(sys.argv[1]+"/"+bin_file):
			if line[0]!=">":
				continue
			contig_bins[line[1:-1]]=bin_name

# store the read names and what bins they belong in in these dictionaries
# strict stores only perfectly aligning reads and permissive stores any aligned reads

print "Parsing sam file and writing reads to appropriate files depending what bin they alligned to..."
files={}
opened_bins={}
for line in sys.stdin:
	if line[0]=="@": 
		continue
	cut = line.strip().split("\t")
	binary_flag = bin(int(cut[1]))

	if binary_flag[-7]=="1":
		F_line=line
		continue
	elif binary_flag[-8]=="1":
		R_line=line

	# get fields for forward and reverse reads	
	F_cut = F_line.strip().split("\t")
	R_cut = R_line.strip().split("\t")

	# skip non aligned reads
	if F_cut[2]=="*" and R_cut[2]=="*":
		continue

	# make sure the R and F reads aligned to the same bin (or allow for one to be unmapped)
	if F_cut[2] != R_cut[2]:
		# check if at least one of the reads is in one of the bins (the other can be unaligned)
		if (F_cut[2] not in contig_bins and F_cut[2]!="*") or (R_cut[2] not in contig_bins and R_cut[2]!="*"):
			continue
		bin1 = contig_bins[F_cut[2]]
		bin2 = contig_bins[R_cut[2]]
		if bin1 != bin2 and bin1 != "unalinged" and bin2 != "unalinged": 
			continue
		bin_name=bin1
	else:
		contig=F_cut[2]
		if contig not in contig_bins:
			continue
		bin_name = contig_bins[contig]

	# open the revelant output files
	if bin_name not in opened_bins:
		opened_bins[bin_name]=None
		files[sys.argv[2]+"/"+bin_name+".R1.fastq.gz"]=gzip.open(sys.argv[2]+"/"+bin_name+".R1.fastq.gz", "wb")
		files[sys.argv[2]+"/"+bin_name+".R2.fastq.gz"]=gzip.open(sys.argv[2]+"/"+bin_name+".R2.fastq.gz", "wb")

	# determine alignment type from bitwise FLAG
	F_binary_flag = bin(int(F_cut[1]))
	R_binary_flag = bin(int(R_cut[1]))

	# if the reads are reversed, fix them
	if F_binary_flag[-5]=='1':
		F_cut[9] = rev_comp(F_cut[9])
		F_cut[10] = F_cut[10][::-1]
	if R_binary_flag[-5]=='1':
		R_cut[9] = rev_comp(R_cut[9])
		R_cut[10] = R_cut[10][::-1]

	# strict assembly
	files[sys.argv[2]+"/"+bin_name+".R1.fastq.gz"].write('@' + F_cut[0] + "/1" + "\n" + F_cut[9] + "\n+\n" + F_cut[10] + "\n")
	files[sys.argv[2]+"/"+bin_name+".R2.fastq.gz"].write('@' + R_cut[0] + "/2" + "\n" + R_cut[9] + "\n+\n" + R_cut[10] + "\n")


print "closing files"
for f in files:
	files[f].close()


print "Finished splitting reads!"


