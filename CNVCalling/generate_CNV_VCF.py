#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2021/10/19
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script transform CNV files obtained with Control FREEC to a gVCF matrix.
It takes as input :
	-a: the anno table used by Control FREEC
		File name must begin with StrainName.xxxx
		StrainName will be the header of each column in the VCF
	-f: the output files .annos.tsv of Control FREEC to be transformed
	-r: the reference genome
	-o: output file
'''
# ---------------------------------------------------------------------------
import csv
import os
import sys
import pandas as pd
import numpy as np
import argparse
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-a", "--anno_table", help="Annotation Table as required for Control FREEC", required=True)
parser.add_argument("-f", "--CNV_files", help="tsv output of Control FREEC containing CNVs", nargs='+', required=True)
parser.add_argument("-r", "--reference", help="fasta file of the reference genome", required=True)
parser.add_argument("-o", "--output", help="Name of the outputed VCF", required=True)

# Read arguments from the command line
args = parser.parse_args()

annoTablePath=args.anno_table
CNVfiles=args.CNV_files
nbFiles=str(len(CNVfiles))
refPath=args.reference
outputPath=args.output
outputTemporary="temporary.vcf"

print('\nRUNNING SCRIPT WITH OPTIONS:')
print('\tAnnotation table:')
print('\t\t'+annoTablePath)
print('\tCNV files ('+nbFiles+' files):')
for element in CNVfiles:
	print('\t\t'+element)
print('\tReference genome:')
print('\t\t'+refPath+'\n')
print('\tOutput:')
print('\t\t'+outputPath+'\n')

# =================================
# Create VCF with all possible CNVs
# =================================

# Read reference to write header
refChr=[]
refSeq=[]
seq=""
ref=open(refPath, 'r')
for line in ref.readlines():
	if line.startswith(">"):
		refChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
		if seq != "":
			refSeq += [seq]
		seq=""
	else:
		seq += line.strip()
refSeq += [seq]
ref.close()
refLen=[len(x) for x in refSeq]

# Create header
VCF=pd.DataFrame([["##fileformat=VCFv4.2"]])
for i in range(len(refChr)):
	VCF=VCF.append([["##contig=<ID=" + refChr[i] + ",length=" + str(refLen[i]) + ">"]])
VCF=VCF.append([["##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"]])
VCF=VCF.append([["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]])

# Read the anno table and creat VCF
tsv_file=open(annoTablePath)
annoTable = csv.reader(tsv_file, delimiter="\t")

for row in annoTable:
	if row[0] != "" and row[4] == "CDS":
		CHROM=row[1].split(".")[0]
		POS=row[5]
		ID=row[3]
		REF="N" # reference: Normal
		ALT="ALT"
		QUAL="."
		FILTER="."
		INFO="."
		FORMAT="GT"

		VCF=VCF.append([[CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]])
tsv_file.close()

# Put index number at 0 1 2 3 4 ...
VCF=VCF.set_axis(list(range(VCF.shape[0])), axis="index")

# ===============================
# Add the genotype of each strain
# ===============================

nbProcessed=0
for CNVfile in CNVfiles:

	# ========================
	# Create the column to add
	# ========================

	# Find the strain name
	strain=CNVfile.split('/')[-1].split('.')[0]

	# Header of the col, that is a blank and the strain name
	head = pd.DataFrame({"CNK": [None]*(1+len(refChr)+1) + [strain]})

	# Get the CNV genotypes
	# Read the CNV file
	values = pd.read_csv(CNVfile, sep="\t")
	# Put variant to "expected" if coverage is lower than 50%
	cond = values["Coverage"] <= 50
	values.loc[cond, "CNK"] = "expected"
	# Get only "CNK" column
	values=pd.DataFrame(values["CNK"])

	# Combine genotype and header and change index labels
	col2add = head.append(values)
	col2add=col2add.set_axis(list(range(col2add.shape[0])), axis="index")
	col2add = col2add.rename(columns={'CNK': strain})

	# Concatenate with the VCF previously created and change col label
	VCF = pd.concat([VCF, col2add], axis=1)
	VCF=VCF.set_axis(list(range(VCF.shape[1])), axis="columns")

	nbProcessed+=1
	print(strain+" processed ("+str(nbProcessed)+"/"+nbFiles+" done)")

# Write to VCF temporary file
VCF.to_csv(outputTemporary, sep="\t", header=False, index=False)

# =======================================
# Split each row in gain and loss variant
# =======================================

print("\nDiscriminating gain and loss variants\n")


# Write to output file
output=open(outputPath, 'w')

# Open the VCF temporary file
tsv_file=open(outputTemporary)
VCF = csv.reader(tsv_file, delimiter="\t")
for row in VCF:
	if row[0].startswith("#"):
		output.write('\t'.join(row)+'\n')
	elif 'loss' not in row and 'gain' not in row:
		pass # If there is no gain or no loss in any strain for a variant, the variant is not written
	else:
		# Write the gain variant
		gainrow=row.copy()
		gainrow[2]="gain_"+gainrow[2]
		gainrow[4]="G" # alternative: Gain
		gainrow = np.array(gainrow)
		gainrow = np.where(gainrow == 'expected', '0/0', gainrow)
		gainrow = np.where(gainrow == 'loss', '0/0', gainrow)
		gainrow = np.where(gainrow == 'gain', '1/1', gainrow)
		# only write if not only 0 genotypes
		if sum(np.char.replace(gainrow[9:], '/', '').astype('int')) > 0 :
			output.write('\t'.join(gainrow)+'\n')

		# Write the loss variant
		lossrow=row.copy()
		lossrow[2]="loss_"+lossrow[2]
		lossrow[4]="L" # alterntive: Loss
		lossrow = np.array(lossrow)
		lossrow = np.where(lossrow == 'expected', '0/0', lossrow)
		lossrow = np.where(lossrow == 'loss', '1/1', lossrow)
		lossrow = np.where(lossrow == 'gain', '0/0', lossrow)
		# only write if not only 0 genotypes
		if sum(np.char.replace(lossrow[9:], '/', '').astype('int')) > 0 :
			output.write('\t'.join(lossrow)+'\n')

output.close()
print("\nVCF file created succesfully")
print("\nFile saved in: "+outputPath)

# Removing VCF temporary file
os.remove(outputTemporary)


