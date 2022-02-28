###################################################################################
## This script was written to polish long-read assemblies with long reads using
## Medaka (https://github.com/nanoporetech/medaka) run on CPUs.
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Run medaka to polish draft assemblies 
conda activate medaka_v1_5
conda list
mkdir 13.medaka

for GENOME in *_sample00_flye.fasta
do 
	HEADER=$(echo ${GENOME} | sed 's/_sample00_flye.fasta//');
	medaka_consensus -i 05.seqtk/${HEADER}_trimfilt5.fastq.gz -d ${GENOME} -o 13.medaka/${HEADER}_medaka -t $THREADS -m r941_min_sup_g507;
	mv ${GENOME} 13.medaka/;
	mv 13.medaka/${HEADER}_medaka/consensus.fasta 13.medaka/${HEADER}_medaka/${HEADER}_flye00_medakaconsensus.fasta;
	mv 13.medaka/${HEADER}_medaka/${HEADER}_flye00_medakaconsensus.fasta .
done

rm *.fai
rm *.mmi
echo "Medaka polishing complete"
conda deactivate
