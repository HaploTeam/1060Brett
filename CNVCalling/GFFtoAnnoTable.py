#!/usr/bin/env python
# -*- coding: utf8 -*-
import sys

# This script convert a GFF annotation file to the tsv file required by CNFREEC
# It takes the GFF file as input



def GFF2AnnoTable(GFFPath):

	outPath=GFFPath.split('.gff')[0]+"_anno_table.tsv"
	outFile=open(outPath, 'w')
	outFile.write("\tContig\tEnd\tID\tKind\tStart\n")

	GFF=open(GFFPath, 'r')
	annotationNb=0
	for line in GFF.readlines():
		if not line.startswith('#'):
			Contig=line.split('\t')[0]
			End=line.split('\t')[4]
			ID=line.split('\n')[0].split('\t')[8].split(';')[0].split('=')[1]+'_'+line.split('\t')[2]
			Kind=line.split('\t')[2]
			Start=line.split('\t')[3]

			outFile.write(str(annotationNb)+'\t'+Contig+'\t'+End+'\t'+ID+'\t'+Kind+'\t'+Start+'\n')

			annotationNb+=1





if __name__ == "__main__":

	GFFPath = sys.argv[1]
	GFF2AnnoTable(GFFPath)

	