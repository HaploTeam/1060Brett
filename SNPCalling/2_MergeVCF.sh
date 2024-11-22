#!/bin/bash
#SBATCH -p fast          # The account name for the job
#SBATCH --job-name=",jobname2,"-All2nOnly-full  # The job name
#SBATCH -o ",scripts_dir,"logs/",jobname2,"-All2nOnly-full.out
#SBATCH -e ",scripts_dir,"logs/",jobname2,"-All2nOnly-full.err
#SBATCH -c 20                 # The number of cpu cores to use
#SBATCH --time=119:59:00       # The time the job will take to run 
#SBATCH --mem=4gb

date
source ~/anaconda3/etc/profile.d/conda.sh
conda activate bwaetc


for i in *.vcf.gz
do echo $i
bcftools merge --threads 20 -Oz -o ${i}.vcf.gz $(cat $i | tr '\\n' ' ')
bcftools index -f ${i}.vcf.gz
done

bcftools merge --threads 20 -Oz -o All2nOnly-full.vcf.gz *vcf.gz
bcftools index -f All2nOnly-full.vcf.gz
bcftools +fill-tags All2nOnly-full.vcf.gz -- -t AF,AN,AC | bgzip > All2nOnly-full.ANAFAC.vcf.gz
bcftools index -f All2nOnly-full.ANAFAC.vcf.gz
bcftools view -i '(GQ>=20 | GQ == \".\") & AC > 0' --types snps -Oz -o All2nOnly-full.ANAFAC.gq20.SNPs.var.vcf.gz All2nOnly-full.ANAFAC.vcf.gz
bcftools index -f All2nOnly-full.ANAFAC.gq20.SNPs.var.vcf.gz
