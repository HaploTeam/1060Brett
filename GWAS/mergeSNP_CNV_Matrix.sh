#!/bin/bash
#MSUB -r genMatrix
#MSUB -n 1
#MSUB -c 1
#MSUB -T 3600
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

module load bcftools

# SNP matrix with all SNPs
PLINK_SNP=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/AllSNPsMatrix/AllSNPs.GQ20.DP5.cov5pto95p.Miss10-Mind20.Biallelic.MAF5.plink
# CNV matrix
PLINK_CNV=/ccc/scratch/cont007/fg0006/loeglerv/BrettGWAS/CNVMatrix/CNVmatrix.1026Strains.MAF5.plink
# Output
PLINK_MERGED=SNP_CNV_matrix.1026Strains.MAF5.plink

# Merge SNPs and CNVs matrices
$CCCSCRATCHDIR/plink --bfile $PLINK_SNP --bmerge $PLINK_CNV --out $PLINK_MERGED --make-bed
