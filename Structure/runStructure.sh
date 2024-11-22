#!/bin/bash
#SBATCH -p fast               
#SBATCH -n 30
#SBATCH --mem-per-cpu=1000
#SBATCH -o slurm.%N.%j.out        
#SBATCH -e slurm.%N.%j.err   

# Convert VCF to plink
plink --vcf $VCF --out $PLINK_OUT.plink --make-bed --allow-extra-chr --const-fid

source /home/vloegler/.bashrc
conda activate FastStructureEnv


# FastStructure command

for K in $(seq 1 15)
do
	structure.py -K $K --input=$WHOLEGENOMESNPs.plink --output=$OUTPUT.FastStruture --seed=100
done

# Admixture command

for K in $(seq 1 15)
do
admixture --cv $PRIMARYSNPs.plink.bed $K | tee log.out
admixture --cv $ACQUIREDSNPs.plink.bed $K | tee log.out
done

