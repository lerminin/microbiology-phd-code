#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=idba
#SBATCH --error=idba.err
#SBATCH --output=idba.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256000M

# Load environments:
module load StdEnv/2020 idba-ud/1.1.3

idba

for R1 in /home/scratch/210610_Illumina/*_R12.fasta
do
        HEADER=$(echo ${R1} | sed 's/_R12.fasta//');
        idba -l $R1 --pre_correction -o ${HEADER}_idba
done

echo "idba-ud complete"
