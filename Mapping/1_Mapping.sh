#!/bin/bash
#SBATCH -p fast          # The account name for the job
#SBATCH --job-name=",jobname,"  # The job name
#SBATCH -o bwa.out
#SBATCH -e bwa.err
#SBATCH -c 5                 # The number of cpu cores to use
#SBATCH --time=11:59:00       # The time the job will take to run 
#SBATCH --mem=5gb

date
source ~/anaconda3/etc/profile.d/conda.sh
conda activate bwaetc

bwa mem -t 5 -M -T 0 -a $REF $READS1 READS2 > $SAMPLE.sam 
samtools sort -m 800M -@ 5 -O BAM -o $SAMPLE.bam $SAMPLE.sam
samtools index $SAMPLE.bam
rm -f $SAMPLE.sam