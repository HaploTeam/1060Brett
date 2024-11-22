#!/bin/bash
#SBATCH -p fast          # The account name for the job
#SBATCH --job-name=",jobname,"  # The job name
#SBATCH -o ",scripts_dir,"logs/",jobname,"-",basename(bams_df$bam[bamrow]),".out
#SBATCH -e ",scripts_dir,"logs/",jobname,"-",basename(bams_df$bam[bamrow]),".err
#SBATCH -c 1                 # The number of cpu cores to use
#SBATCH --time=11:59:00       # The time the job will take to run 
#SBATCH --mem=4gb

date
source ~/anaconda3/etc/profile.d/conda.sh
conda activate bwaetc

bcftools mpileup --per-sample-mF --redo-BAQ --min-BQ 30 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --threads 1 -Ob -o $SAMPLE.mpileup -f $Reference.fasta $SAMPLE.bam
bcftools call --threads 1 -a GQ --ploidy $PLOIDY -m -Oz -o $SAMPLE.mpileup.vcf.gz $SAMPLE.mpileup
bcftools index -f $SAMPLE.mpileup.vcf.gz && rm -f $SAMPLE.mpileup
bcftools view -i 'QUAL>99' -Oz -o $SAMPLE.mpileup.QUAL100.vcf.gz $SAMPLE.mpileup.vcf.gz
bcftools index -f $SAMPLE.mpileup.QUAL100.vcf.gz
bgzip -t $SAMPLE.mpileup.vcf.gz 2> $SAMPLE.mpileup.vcf.gz.t
rm -f $SAMPLE.mpileup.vcf.gz