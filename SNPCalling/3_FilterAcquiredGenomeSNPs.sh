#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1
#SBATCH -c 4                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=400              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load bcftools

WORKDIR=/shared/home/vloegler/BrettDiversity/04-analysis
DATADIR=/shared/home/vloegler/BrettDiversity/03-data
SCRIPTDIR=/shared/home/vloegler/BrettDiversity/02-scripts

mkdir $WORKDIR/03_AcquiredGenomeSNPs
cd $WORKDIR/03_AcquiredGenomeSNPs

# Filter Acquired genome SNPs

# 1. Change sample names
# 2. Filter 50% non missing
# 3. SubSet non-admixed samples

# Get samples new names
bcftools query -l $DATADIR/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.vcf.gz | cut -d '/' -f 7 | cut -d '-' -f 1 > NewSamples.txt

# Get list of non admixed samples
cat $DATADIR/Clusters.txt | cut -f 1 > NonAdmixedSamples.txt

bcftools reheader --sample NewSamples.txt $DATADIR/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.vcf.gz | \
	bcftools view -i 'F_MISSING<0.5' | \
	bcftools view --samples-file NonAdmixedSamples.txt --force-samples | \
	bcftools +fill-tags | \
	bcftools view --threads 4 -Oz -o Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.Miss50-Mind70.NonAdmixed.vcf.gz

