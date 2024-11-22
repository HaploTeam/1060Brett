import sys
import os
import glob

import pandas as pd
pd.set_option('display.width', 1000)

workingPath = os.getcwd()
CNFreecPath = "/".join(workingPath.split('/')[:-1])

sys.path.insert(0, CNFreecPath + "/MappingTools-master/mappingTools")
# sys.path.insert(0, "/ccc/work/cont007/fg0006/friedria/Scripts")
sys.path.insert(0, "/home/vloegler/Scripts")

import files
from wrappers.freec import launch_bam
from wrappers.processing import multiprocess_this, subprocess_this
from freec_analysis import process_fname

from collections import Counter

dname = os.path.dirname
bname = os.path.basename

freec_p = "/home/vloegler/Scripts/FREEC-11.6/src/freec"
samtools_p = "samtools"
sample_funname = lambda fname : bname(fname).split(".")[0]

def get_file(name) :
	fname = os.path.realpath(__file__)
	return os.path.join(dname(fname), name)

def launch_freec_args(fname, freec_p, samtools_p, reference, outdir, ignore_exist, config, annos, defaultcn) :

	# launch freec
	launch_bam(fname, freec_p, samtools_p, reference, ignore_exist=ignore_exist, config=config, outdir=outdir)

	# analyse result
	outfile = os.path.join(outdir, bname(fname) + "_CNVs")
	if not os.path.isfile(outfile) :
		print ("Fail for : %s" %(fname))
		return None

	return process_fname(outfile, annos, defaultcn, ignore_exist=False)

def launch(bamfiles, ncore=20, ignore_exist=True) :

	reference = get_file("Bbruxellensis.fasta")
	mingc = 0.3
	maxgc = 0.5

	ploidy = get_file("ploidy.csv")
	ploidy = pd.read_csv(ploidy, sep=";", index_col=0)
	ploidy = ploidy.set_index("Standardized name")["Ploidy"].to_dict() 
	ploidy = {sample : int(value) for sample, value in ploidy.items()}

	bamfiles = glob.glob(bamfiles) if isinstance(bamfiles, str) else bamfiles
	args = []

	annos = get_file("Brett_anno_table.tsv")
	annos = pd.read_csv(annos, sep="\t", index_col=0)
	
	# we remove the gff contig feature and remove the chr part due to Freec BUG
	# WARNING : IN THIS SCRIPT WE ONLY LOOK AT CDS. CHANGE IF NEED

	annos = annos[annos["Kind"] == "CDS"]
	#annos["Contig"] = annos["Contig"].str.slice(start=3) # in comment because here contigs are not named 'chrX'
	annos["Contig"] = annos["Contig"].str.slice(stop=7)
	annos = annos[annos["End"] - annos["Start"] > 0]

	for fname in bamfiles :

		sample = sample_funname(fname)
		sample_ploidy = ploidy.get(sample)
		if not sample_ploidy : print ("Ploidy not found for : %s (%s)" %(sample, fname)); continue

		config = {"sample" : {"mateOrientation" : "FR"}, "general" : {"ploidy" : sample_ploidy,
		"minExpectedGC" : "%.1f" %(mingc), "maxExpectedGC" : "%.1f" %(maxgc)}}
	
		soutdir = fname + ".freec"
		args.append((fname, freec_p, samtools_p, reference, soutdir, ignore_exist, config, annos, sample_ploidy))

	args = sorted(args, key = lambda x : x[0])
	print(args)
	res = multiprocess_this(launch_freec_args, args, ncore=ncore)


args = sys.argv[1:]
all_args = files.get_args(args)
try:
    bam = all_args["b"]
except:
    bam = ""

if __name__ == "__main__" :
	resfiles = launch(bam, ncore=1)
