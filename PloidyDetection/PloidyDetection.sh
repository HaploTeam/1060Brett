#!/bin/bash
#SBATCH -p fast          # The account name for the job
#SBATCH --job-name=",jobname,"  # The job name
#SBATCH -o ",scripts_dir,"logs/",jobname,"-$SAMPLE.nquire.out
#SBATCH -e ",scripts_dir,"logs/",jobname,"-$SAMPLE.nquire.err
#SBATCH -c 1                 # The number of cpu cores to use
#SBATCH --time=11:59:00       # The time the job will take to run 
#SBATCH --mem=4gb

date
cd /shared/home/jnrunge/data/bruxellensis/aggregated/AD_nQuire/
~/software/nQuire/nQuire/nQuire create -c 6 -b $SAMPLE -o $SAMPLE.nquire
~/software/nQuire/nQuire/nQuire create -c 6 -b $SAMPLE -o $SAMPLE.nquire-nocutoff -f 0
~/software/nQuire/nQuire/nQuire denoise $SAMPLE.nquire.bin -o $SAMPLE.nquire-denoised
rm -f $SAMPLE.nquire.bin
~/software/nQuire/nQuire/nQuire histotest $SAMPLE.nquire-denoised.bin > $SAMPLE.nquire-denoised-histo.test
