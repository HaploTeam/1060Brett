#!/bin/bash
#SBATCH -p fast          # The account name for the job
#SBATCH --job-name=",jobname,"  # The job name
#SBATCH -o ",scripts_dir,"logs/",basename(bam),"-extract1n2n.out
#SBATCH -e ",scripts_dir,"logs/",basename(bam),"-extract1n2n.err
#SBATCH -c 1                 # The number of cpu cores to use
#SBATCH --time=11:59:00       # The time the job will take to run 
#SBATCH --mem=1gb

source ~/anaconda3/etc/profile.d/conda.sh
conda activate bwaetc

samtools view -O BAM -o $SAMPLE.2n.mq30.bam --region-file Reference.2ncontigs -q 30 $SAMPLE.bam
samtools view -O BAM -o $SAMPLE.1n.mq30.bam --region-file Reference.1ncontigs -q 30 $SAMPLE.bam
