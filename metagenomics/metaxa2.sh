#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=metaxa2
#SBATCH --error=metaxa2.err
#SBATCH --output=metaxa2.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256000M

# Load environments:
module load StdEnv/2020  gcc/9.3.0 metaxa2/2.2

metaxa2 -h

for R1 in /home/scratch/210610_Illumina/*R1_illumina_cleantrim.fastq.gz
do
        HEADER=$(echo ${R1} | sed 's/_R1_illumina_cleantrim.fastq.gz//');
        R2=$(echo ${R1} | sed 's/R1/R2/')
        metaxa2 -1 -2 -o ${HEADER}_metaxa2 --cpu 32 -f q -g ssu
        mkdir ${HEADER}_metaxa2_ssu
        mv ${HEADER}_metaxa2.* ${HEADER}_metaxa2_ssu
done

for R1 in /home/scratch/210610_Illumina/*R1_illumina_cleantrim.fastq.gz
do
        HEADER=$(echo ${R1} | sed 's/_R1_illumina_cleantrim.fastq.gz//');
        R2=$(echo ${R1} | sed 's/R1/R2/')
        metaxa2 -1 -2 -o ${HEADER}_metaxa2 --cpu 32 -f q -g lsu
        mkdir ${HEADER}_metaxa2_lsu
        mv ${HEADER}_metaxa2.* ${HEADER}_metaxa2_lsu
done

echo "metaxa2 complete"
