#!/bin/bash
#SBATCH -p fast               # partition all / fast
#SBATCH -n 40                      # number of cores
#SBATCH -o slurm.%N.%j.out        # STDOUT
#SBATCH -e slurm.%N.%j.err        # STDERR

function run_parallel ()
{
srun -n1 --exclusive "$@"
}

module load samtools 

#for bam in /home/vloegler/BAMfiles_Brett/YJS8833.rg.dedup.bam
for bam in /home/vloegler/BAMfiles_Brett/*.bam
do
    run_parallel python3 cnfreec_brett.py -b $bam &

done
wait