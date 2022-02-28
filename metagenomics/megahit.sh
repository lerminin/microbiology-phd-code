#!/bin/bash

#SBATCH --time=30:00:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=megahit
#SBATCH --error=megahit.err
#SBATCH --output=megahit.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000M

# Load environments:
module load StdEnv/2020 megahit

megahit -h

for R1 in /home/scratch/210610_Illumina/*R1_illumina_cleantrim.fastq.gz
do
	HEADER=$(echo ${R1} | sed 's/_R1_illumina_cleantrim.fastq.gz//');
	R2=$(echo ${R1} | sed 's/R1/R2/')
	megahit -1 $R1 -2 $R2 -t 32 --presets meta-large -o ${HEADER}_megahit
done

echo "megahit complete" 
