#!/bin/bash
#MSUB -r filterVCF
#MSUB -n 1
#MSUB -c 4
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,store,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools

FULLMATRIX=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/AllSNPsMatrix/All2nOnly-All2nOnly-full.ANAFAC.SNPs.var.GQ20.DP5.cov5pto95p.annotatedID.vcf.gz
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/AllSNPsMatrix/
PREFIX=All2nOnly-All2nOnly-full.ANAFAC.SNPs.var.GQ20.DP5.cov5pto95p.annotatedID

# Rename and filter for Miss10
bcftools query -l $FULLMATRIX | cut -d "/" -f 7 | cut -d "-" -f 1 > NewNames.txt
bcftools reheader --samples NewNames.txt $FULLMATRIX | bcftools view -i "F_MISSING<0.1" -Oz -o $PREFIX.Miss10.vcf.gz --threads 8
rm -f NewNames.txt

# Filter Mind20
$CCCSCRATCHDIR/plink --vcf $PREFIX.Miss10.vcf.gz --mind 0.2 --const-fid --allow-extra-chr --make-bed --out $PREFIX.Miss10-Mind20-plink
cut -f 2 $PREFIX.Miss10-Mind20-plink.irem > $PREFIX.Miss10-Mind20-plink.irem.samples
bcftools view -S ^$PREFIX.Miss10-Mind20-plink.irem.samples -Oz -o $PREFIX.Miss10-Mind20.vcf.gz $PREFIX.Miss10.vcf.gz
rm -f $PREFIX.Miss10-Mind20-plink.*

# Filter biallelic, filter MAF, add rs ID and change chromosome name
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
bcftools view -m2 -M2 $PREFIX.Miss10-Mind20.vcf.gz | bcftools view --min-af 0.05:minor  | bcftools annotate --set-id 'rs_%CHROM\_%POS' | bcftools annotate --rename-chrs $CHRNAMES --threads 4 -Oz -o $PREFIX.Miss10-Mind20.Biallelic.MAF5.plinkReady.vcf.gz
rm -f $CHRNAMES

# Convert to plink
$CCCSCRATCHDIR/plink --vcf $PREFIX.Miss10-Mind20.Biallelic.MAF5.plinkReady.vcf.gz --out AllSNPs.GQ20.DP5.cov5pto95p.Miss10-Mind20.Biallelic.MAF5.plink --make-bed
