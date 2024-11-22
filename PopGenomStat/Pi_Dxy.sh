#!/bin/bash
#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 3
#SBATCH -c 2                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=400              # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o slurm.%N.%j.out      # File to which STDOUT will be written
#SBATCH -e slurm.%N.%j.err      # File to which STDERR will be written

module load bcftools

source /shared/home/vloegler/.bashrc
conda activate PixyEnv

WORKDIR=/shared/home/vloegler/BrettDiversity/04-analysis
DATADIR=/shared/home/vloegler/BrettDiversity/03-data
SCRIPTDIR=/shared/home/vloegler/BrettDiversity/02-scripts

#mkdir $WORKDIR/04_PixyResults
cd $WORKDIR/04_PixyResults

# Index vcf files
#srun --exclusive -N1 -n1 bcftools index -t $WORKDIR/01_WholeGenomeSNPs/All2nOnly-All2nOnly-full.ANAFAC.SNPs.Novar.GQ20.DP5.cov5pto95p.Miss10-Mind20.NonAdmixed.vcf.gz &
#srun --exclusive -N1 -n1 bcftools index -t $WORKDIR/02_PrimaryGenomeSNPs/All2n-v5.ANAFAC.SNPs.Novar.GQ20.DP5.cov5pto95p.Miss10-Mind20.NonAdmixed.vcf.gz &
#srun --exclusive -N1 -n1 bcftools index -t $WORKDIR/03_AcquiredGenomeSNPs/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.Miss50-Mind70.NonAdmixed.vcf.gz &
#wait

# create population file for allopolyploids
bcftools query -l $WORKDIR/03_AcquiredGenomeSNPs/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.Miss50-Mind70.NonAdmixed.vcf.gz > AllopolyploidIsolates.txt
grep -f AllopolyploidIsolates.txt $DATADIR/Clusters.txt > ClustersAllopolyploid.txt

# Compute pi and dxy per clade
#srun --exclusive -N1 -n1 pixy --stats pi dxy \
#							--vcf $WORKDIR/01_WholeGenomeSNPs/All2nOnly-All2nOnly-full.ANAFAC.SNPs.Novar.GQ20.DP5.cov5pto95p.Miss10-Mind20.NonAdmixed.vcf.gz \
#							--populations $DATADIR/Clusters.txt \
#							--window_size 10000 \
#							--output_prefix Whole \
#							--n_cores 4 &
srun --exclusive -N1 -n1 pixy --stats pi dxy \
							--vcf $WORKDIR/02_PrimaryGenomeSNPs/All2n-v5.ANAFAC.SNPs.Novar.GQ20.DP5.cov5pto95p.Miss10-Mind20.NonAdmixed.vcf.gz \
							--populations $DATADIR/Clusters.txt \
							--window_size 10000 \
							--output_prefix Primary \
							--n_cores 4 &
#srun --exclusive -N1 -n1 pixy --stats pi dxy \
#							--vcf $WORKDIR/03_AcquiredGenomeSNPs/Merge-extracted1ns.ANAFAC.SNPs.GQ20.DP3.cov5pto95p.Miss50-Mind70.NonAdmixed.vcf.gz \
#							--populations ClustersAllopolyploid.txt \
#							--window_size 10000 \
#							--output_prefix Acquired \
#							--n_cores 4 &
wait

sed 's/^/Whole\t/g' Whole_pi.txt | sed '0,/Whole/{s/Whole/Genome/}' > tmp1
sed 's/^/Primary\t/g' Primary_pi.txt | sed '0,/Primary/{s/Primary/Genome/}' | grep -v Genome > tmp2
sed 's/^/Acquired\t/g' Acquired_pi.txt | sed '0,/Acquired/{s/Acquired/Genome/}' | grep -v Genome > tmp3
cat tmp1 tmp2 tmp3 > Brbr.946Isolates.WholePrimaryAcquired.Pixy.Pi.tsv
sed 's/^/Whole\t/g' Whole_dxy.txt | sed '0,/Whole/{s/Whole/Genome/}' > tmp1
sed 's/^/Primary\t/g' Primary_dxy.txt | sed '0,/Primary/{s/Primary/Genome/}' | grep -v Genome > tmp2
sed 's/^/Acquired\t/g' Acquired_dxy.txt | sed '0,/Acquired/{s/Acquired/Genome/}' | grep -v Genome > tmp3
cat tmp1 tmp2 tmp3 > Brbr.946Isolates.WholePrimaryAcquired.Pixy.Dxy.tsv


