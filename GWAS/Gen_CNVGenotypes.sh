#!/bin/bash
#MSUB -r SNPfiltering
#MSUB -n 1
#MSUB -c 1
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools

# Filtering CNV for AllSNPs GWAS
PREFIX=CNV_matrix_total
ALLSNPSMATRIX=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/AllSNPsMatrix/All2nOnly-All2nOnly-full.ANAFAC.SNPs.var.GQ20.DP5.cov5pto95p.annotatedID.Miss10-Mind20.Biallelic.MAF5.plinkReady.vcf.gz

# Get 1026 strains
bcftools query -l $ALLSNPSMATRIX > 1026Strains.txt

# change chr names, filter samples and MAF 5%
CHRNAMES=ChrNames.txt
rm -f $CHRNAMES
echo "DEBR0S1 1" >> $CHRNAMES
echo "DEBR0S2 2" >> $CHRNAMES
echo "DEBR0S3 3" >> $CHRNAMES
echo "DEBR0S4 4" >> $CHRNAMES
echo "DEBR0S5 5" >> $CHRNAMES
echo "DEBR0S6 6" >> $CHRNAMES
echo "DEBR0S7 7" >> $CHRNAMES
echo "DEBR0S8 8" >> $CHRNAMES
bcftools annotate --rename-chrs $CHRNAMES $PREFIX.vcf | bcftools view -S 1026Strains.txt | bcftools view --min-af 0.05:minor -o $PREFIX.1026Strains.MAF5.plinkReady.vcf
rm -f $CHRNAMES

# Convert to plink
$CCCSCRATCHDIR/plink --vcf $PREFIX.1026Strains.MAF5.plinkReady.vcf --out CNVmatrix.1026Strains.MAF5.plink --make-bed
