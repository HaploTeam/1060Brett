#!/bin/bash
#SBATCH -p fast          # The account name for the job
#SBATCH --job-name=",jobname,"  # The job name
#SBATCH -o ",scripts_dir,"logs/",basename(bam),"-",jobname,".out
#SBATCH -e ",scripts_dir,"logs/",basename(bam),"-",jobname,".err
#SBATCH -c 2                 # The number of cpu cores to use
#SBATCH --time=11:59:00       # The time the job will take to run 
#SBATCH --mem=3gb

conda activate bwaetc

samtools sort -m 1G -@ 2 -n -O bam -o $SAMPLE.sortn.bam $SAMPLE.bam

samtools fixmate -@ 2 -m $SAMPLE.sortn.bam $SAMPLE.sortn.m.bam

rm -f $SAMPLE.sortn.bam

samtools sort -m 1G -@ 2 -O bam -o $SAMPLE.sortn.m.sortc.bam $SAMPLE.sortn.m.bam

rm -f $SAMPLE.sortn.m.bam

samtools markdup -l $rl -s -@ 2 $SAMPLE.sortn.m.sortc.bam $SAMPLE.sortn.m.sortc.markdup.bam

rm -f $SAMPLE.sortn.m.sortc.bam

samtools index $SAMPLE.sortn.m.sortc.markdup.bam