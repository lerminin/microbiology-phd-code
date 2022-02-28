#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=plasspades
#SBATCH --error=pspades.err
#SBATCH --output=pspades.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256000M

# Load environments:
module load StdEnv/2020 spades/3.15.1

spades.py -h

for R1 in /home/scratch/210610_Illumina/*R1_illumina_cleantrim.fastq.gz
do
	HEADER=$(echo ${R1} | sed 's/_R1_illumina_cleantrim.fastq.gz//');
	R2=$(echo ${R1} | sed 's/R1/R2/')
	spades.py -1 $R1 -2 $R2 -t 32 --meta -o ${HEADER}_pspades --plasmid
done

echo "spades complete" 
