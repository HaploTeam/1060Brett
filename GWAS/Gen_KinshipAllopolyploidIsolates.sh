#!/bin/bash
#MSUB -r filterVCF
#MSUB -n 1
#MSUB -c 4
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,store
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools

FULLMATRIX=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/AllSNPsMatrix/All2nOnly-All2nOnly-full.ANAFAC.SNPs.var.GQ20.DP5.cov5pto95p.annotatedID.vcf.gz
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/AllSNPsMatrix/
PREFIX=All2nOnly-All2nOnly-full.ANAFAC.SNPs.var.GQ20.DP5.cov5pto95p.annotatedID

# Rename and filter samples (keep only strains in common from 1n and 2n allopolyploid SNP matrices)
bcftools query -l $FULLMATRIX | cut -d "/" -f 7 | cut -d "-" -f 1 > NewNames.txt
MATRIX1N=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/1nSNPsMatrix/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP1.cov5pto95p.var.Miss10-Mind20.vcf.gz
MATRIX2N=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/2nSNPsMatrix/All2n-v5.ANAFAC.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.vcf.gz
bcftools query -l $MATRIX2N > Strains_2nSNPs_Allopolyploid.txt
bcftools query -l $MATRIX1N > Strains_1nSNPs_Allopolyploid.txt
cat Strains_2nSNPs_Allopolyploid.txt Strains_1nSNPs_Allopolyploid.txt | sort | uniq -d > Strains_Allopolyploid_GWAS.txt
NBstrains=$(cat Strains_Allopolyploid_GWAS.txt | wc -l)
rm -f Strains_2nSNPs_Allopolyploid.txt
rm -f Strains_1nSNPs_Allopolyploid.txt
bcftools reheader --samples NewNames.txt $FULLMATRIX | bcftools view --samples-file Strains_Allopolyploid_GWAS.txt --force-samples | bcftools +fill-tags | bcftools view --min-ac 1 --threads 8 -Oz -o $PREFIX.${NBstrains}Strains.vcf.gz
rm -f NewNames.txt

# Change chr names
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

# Filter Miss10, biallelic, add rs ID, change Chr Names and filter MAF
bcftools view -i "F_MISSING<0.1" $PREFIX.${NBstrains}Strains.vcf.gz | bcftools view -m2 -M2 | bcftools view --min-af 0.05:minor | bcftools annotate --set-id 'rs_%CHROM\_%POS' | bcftools annotate --rename-chrs $CHRNAMES --threads 8 -O z -o $PREFIX.${NBstrains}Strains.Miss10.Biallelic.MAF5.plinkReady.vcf.gz

$CCCSCRATCHDIR/plink --vcf $PREFIX.${NBstrains}Strains.Miss10.Biallelic.MAF5.plinkReady.vcf.gz --out AllSNPs.GQ20.DP5.cov5pto95p.annotatedID.Miss10.2n1nStrains.Biallelic.plink --make-bed
