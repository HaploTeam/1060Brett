#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=400              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load bcftools

WORKDIR=/shared/home/vloegler/BrettDiversity/04-analysis
DATADIR=/shared/home/vloegler/BrettDiversity/03-data
SCRIPTDIR=/shared/home/vloegler/BrettDiversity/02-scripts

mkdir $WORKDIR/07_LD
cd $WORKDIR/07_LD

# Get Beer, Wine and Teq isolates
grep Beer $DATADIR/Clusters.txt | cut -f 1 > BeerIsolates.txt
paste BeerIsolates.txt BeerIsolates.txt > tmp
mv -f tmp BeerIsolates.txt
grep Wine1 $DATADIR/Clusters.txt | cut -f 1 > Wine1Isolates.txt
paste Wine1Isolates.txt Wine1Isolates.txt > tmp
mv -f tmp Wine1Isolates.txt
grep Teq $DATADIR/Clusters.txt | cut -f 1 > TeqIsolates.txt
paste TeqIsolates.txt TeqIsolates.txt > tmp
mv -f tmp TeqIsolates.txt

# Convert to plink
plink --vcf $WORKDIR/02_PrimaryGenomeSNPs/All2n-v5.ANAFAC.SNPs.Novar.GQ20.DP5.cov5pto95p.Miss10-Mind20.NonAdmixed.vcf.gz --min-ac 1 --make-bed --out PrimarySNPs.plink -aec
plink --vcf $WORKDIR/03_AcquiredGenomeSNPs/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.Miss50-Mind70.NonAdmixed.vcf.gz --min-ac 1 --make-bed --out AcquiredSNPs.plink -aec

# Filter for pure allopolyploid clusters
# Keep only regions shared between subgenomes
plink --bfile PrimarySNPs.plink --keep BeerIsolates.txt --extract range $DATADIR/SharedCoordinatesBeer.bed --min-ac 1 --out PrimarySNPs.Beer.SharedPositions.plink --make-bed -aec
plink --bfile AcquiredSNPs.plink --keep BeerIsolates.txt --extract range $DATADIR/SharedCoordinatesBeer.bed --min-ac 1 --out AcquiredSNPs.Beer.SharedPositions.plink --make-bed -aec

plink --bfile PrimarySNPs.plink --keep Wine1Isolates.txt --extract range $DATADIR/SharedCoordinatesWine1.bed --min-ac 1 --out PrimarySNPs.Wine1.SharedPositions.plink --make-bed -aec
plink --bfile AcquiredSNPs.plink --keep Wine1Isolates.txt --extract range $DATADIR/SharedCoordinatesWine1.bed --min-ac 1 --out AcquiredSNPs.Wine1.SharedPositions.plink --make-bed -aec

plink --bfile PrimarySNPs.plink --keep TeqIsolates.txt --extract range $DATADIR/SharedCoordinatesTeq.bed --min-ac 1 --out PrimarySNPs.Teq.SharedPositions.plink --make-bed -aec
plink --bfile AcquiredSNPs.plink --keep TeqIsolates.txt --extract range $DATADIR/SharedCoordinatesTeq.bed --min-ac 1 --out AcquiredSNPs.Teq.SharedPositions.plink --make-bed -aec

# Compute LD
plink --bfile PrimarySNPs.Beer.SharedPositions.plink --r2 gz --ld-window-r2 0 --ld-window-kb 5000 -aec --out PrimarySNPs.Beer.SharedPositions.plink
plink --bfile AcquiredSNPs.Beer.SharedPositions.plink --r2 gz --ld-window-r2 0 --ld-window-kb 5000 -aec --out AcquiredSNPs.Beer.SharedPositions.plink

plink --bfile PrimarySNPs.Wine1.SharedPositions.plink --r2 gz --ld-window-r2 0 --ld-window-kb 5000 -aec --out PrimarySNPs.Wine1.SharedPositions.plink
plink --bfile AcquiredSNPs.Wine1.SharedPositions.plink --r2 gz --ld-window-r2 0 --ld-window-kb 5000 -aec --out AcquiredSNPs.Wine1.SharedPositions.plink

plink --bfile PrimarySNPs.Teq.SharedPositions.plink --r2 gz --ld-window-r2 0 --ld-window-kb 5000 -aec --out PrimarySNPs.Teq.SharedPositions.plink
plink --bfile AcquiredSNPs.Teq.SharedPositions.plink --r2 gz --ld-window-r2 0 --ld-window-kb 5000 -aec --out AcquiredSNPs.Teq.SharedPositions.plink


