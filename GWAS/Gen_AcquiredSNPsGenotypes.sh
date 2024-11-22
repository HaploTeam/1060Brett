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

FULLMATRIX=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/1nSNPsMatrix/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP1.cov5pto95p.var.vcf.gz
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/1nSNPsMatrix/
PREFIX=Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP1.cov5pto95p.var

# Rename and filter for Miss10
bcftools query -l $FULLMATRIX | cut -d "/" -f 7 | cut -d "-" -f 1 > NewNames.txt
bcftools reheader --samples NewNames.txt $FULLMATRIX | bcftools view -i "F_MISSING<0.1" -Oz -o $PREFIX.Miss10.vcf.gz --threads 8
rm -f NewNames.txt

# Filter Mind20
$CCCSCRATCHDIR/plink --vcf $PREFIX.Miss10.vcf.gz --mind 0.2 --const-fid --allow-extra-chr --make-bed --out $PREFIX.Miss10-Mind20-plink
cut -f 2 $PREFIX.Miss10-Mind20-plink.irem > $PREFIX.Miss10-Mind20-plink.irem.samples
bcftools view -S ^$PREFIX.Miss10-Mind20-plink.irem.samples -Oz -o $PREFIX.Miss10-Mind20.vcf.gz $PREFIX.Miss10.vcf.gz
rm -f $PREFIX.Miss10-Mind20-plink.*

# Find common Samples between 1n SNPs and 2n SNPs (Samples set changes because of the Mind20 filter)
# For this step, the matrix 2n with Allopolyploid strains and Miss10-Mind20 must have been generated
MATRIX2N=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/2nSNPsMatrix/All2n-v5.ANAFAC.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.vcf.gz
bcftools query -l $MATRIX2N > Strains_2nSNPs_Allopolyploid.txt
bcftools query -l $PREFIX.Miss10-Mind20.vcf.gz > Strains_1nSNPs_Allopolyploid.txt
cat Strains_2nSNPs_Allopolyploid.txt Strains_1nSNPs_Allopolyploid.txt | sort | uniq -d > Strains_Allopolyploid_GWAS.txt
NBstrains=$(cat Strains_Allopolyploid_GWAS.txt | wc -l)
rm -f Strains_2nSNPs_Allopolyploid.txt
rm -f Strains_1nSNPs_Allopolyploid.txt

# Filter Strains, biallelic, filter MAF, add rs ID and change chromosome name
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
bcftools view -S Strains_Allopolyploid_GWAS.txt $PREFIX.Miss10-Mind20.vcf.gz | bcftools +fill-tags | bcftools view -m2 -M2 | bcftools view --min-af 0.05:minor | bcftools annotate --set-id 'rs_%CHROM\_%POS' | bcftools annotate --rename-chrs $CHRNAMES --threads 4 -Oz -o $PREFIX.Miss10-Mind20.${NBstrains}Strains.Biallelic.MAF5.plinkReady.vcf.gz
rm -f $CHRNAMES

# Convert to plink
$CCCSCRATCHDIR/plink --vcf $PREFIX.Miss10-Mind20.${NBstrains}Strains.Biallelic.MAF5.plinkReady.vcf.gz --out 1nSNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.Biallelic.MAF5.plink --make-bed
