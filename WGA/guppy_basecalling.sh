#!/bin/bash

#SBATCH --job-name=basecalling
#SBATCH --gres=gpu:v100l:2 
#SBATCH --cpus-per-task=8  
#SBATCH --mem=10G  
#SBATCH --time=0-12:00 
#SBATCH --mail-type=ALL
#SBTACH --mail-user=email@email.email
#SBATCH --output=guppy_sup.log
#SBATCH --error=guppy_sup.err

module load StdEnv/2020  gcc/9.3.0  cuda/11.4 ont-guppy/5.0.16
echo "modules loaded"

nvidia-smi
echo "gpu detected"

guppy_basecaller \
--input_path /home/scratch/SUP/pass \
--save_path /home/scratch/SUP/fastq \
--config dna_r9.4.1_450bps_sup.cfg -x cuda:all \
--recursive --flowcell FLO-MIN106 --kit SQK-LSK109 

guppy_barcoder \
--input_path /home/scratch/SUP/fastq/ \
--save_path /home/scratch/SUP/demultiplex_SUP/ \
--barcode_kits EXP-NBD103 -x cuda:all --recursive
echo "guppy finished"

