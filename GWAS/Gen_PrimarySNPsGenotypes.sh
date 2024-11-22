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

FULLMATRIX=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/2nSNPsMatrix/All2n-v5.ANAFAC.vcf.gz
WORKDIR=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/2nSNPsMatrix/
PREFIX=All2n-v5.ANAFAC
MATRIX1N=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/1nSNPsMatrix/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP1.cov5pto95p.var.Miss10.vcf.gz
COV5TO95REGIONS=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/RegionsCov5pto95p/2nSNPs_All2n-coverage-filtered-regions.tsv

# Rename samples, Filter Allopolyploid isolates, variable positions, SNPs and GQ20
bcftools query -l $FULLMATRIX | cut -d "/" -f 7 | cut -d "-" -f 1 > NewNames.txt
bcftools query -l $MATRIX1N > 2n1nStrains.txt
bcftools reheader --samples NewNames.txt $FULLMATRIX | bcftools view -S 2n1nStrains.txt | bcftools +fill-tags | bcftools view --min-ac 1 --threads 8 -O z -o $PREFIX.2n1nStrains.var.vcf.gz
bcftools +fill-tags $PREFIX.2n1nStrains.var.vcf.gz | bcftools view -V indels,mnps,ref,bnd,other --threads 8 -O z -o $PREFIX.2n1nStrains.var.SNPs.vcf.gz
bcftools +setGT $PREFIX.2n1nStrains.var.SNPs.vcf.gz -- -t q -n . -i 'FMT/GQ < 20' | bcftools +fill-tags | bcftools view --min-ac 1 --threads 8 -Oz -o $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.vcf.gz
rm -f NewNames.txt
# Filter Cov5pto95p and Miss10
bcftools index $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.vcf.gz
bcftools view -R $COV5TO95REGIONS $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.vcf.gz | bcftools view -i "F_MISSING<0.1" --threads 8 -Oz -o $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10.vcf.gz

# Filter Mind20
$CCCSCRATCHDIR/plink --vcf $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10.vcf.gz --mind 0.2 --const-fid --allow-extra-chr --make-bed --out $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20-plink
cut -f 2 $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20-plink.irem > $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20-plink.irem.samples
bcftools view -S ^$PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20-plink.irem.samples --threads 8 -Oz -o $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.vcf.gz $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10.vcf.gz
rm -f $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20-plink.*

# Find common Samples between 1n SNPs and 2n SNPs (Samples set changes because of the Mind20 filter)
# For this step, the matrix 1n with Allopolyploid strains and Miss10-Mind20 must have been generated
MATRIX1N=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/1nSNPsMatrix/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP1.cov5pto95p.var.Miss10-Mind20.vcf.gz
bcftools query -l $MATRIX1N > Strains_1nSNPs_Allopolyploid.txt
bcftools query -l $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.vcf.gz > Strains_2nSNPs_Allopolyploid.txt
cat Strains_2nSNPs_Allopolyploid.txt Strains_1nSNPs_Allopolyploid.txt | sort | uniq -d > Strains_Allopolyploid_GWAS.txt
NBstrains=$(cat Strains_Allopolyploid_GWAS.txt | wc -l)
rm -f Strains_2nSNPs_Allopolyploid.txt
rm -f Strains_1nSNPs_Allopolyploid.txt

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
bcftools view -S Strains_Allopolyploid_GWAS.txt $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.vcf.gz | bcftools +fill-tags | bcftools view -m2 -M2 | bcftools view --min-af 0.05:minor | bcftools annotate --set-id 'rs_%CHROM\_%POS' | bcftools annotate --rename-chrs $CHRNAMES --threads 8 -Oz -o $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.${NBstrains}Strains.Biallelic.MAF5.plinkReady.vcf.gz
rm -f $CHRNAMES

# Convert to plink
$CCCSCRATCHDIR/plink --vcf $PREFIX.2n1nStrains.var.SNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.${NBstrains}Strains.Biallelic.MAF5.plinkReady.vcf.gz --out 2nSNPs.GQ20.DP1.cov5pto95p.Miss10-Mind20.Biallelic.MAF5.plink --make-bed
