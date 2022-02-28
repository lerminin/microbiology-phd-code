###################################################################################
## This script polishes long read assemblies with Racon (https://github.com/lbcb-sci/racon),
## which is typically run after Medaka polishing. 
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=12

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Generate initial .sam file
conda activate minimap2_samtools

for GENOME in *_trimfilt5.fastq.gz
do 
    HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
    minimap2 -ax map-ont -t $THREADS ${HEADER}_assembly.fasta ${GENOME} > ${HEADER}_assemblymapped.sam
done
conda deactivate

# Run Racon for the 1st time
conda activate racon
conda list

for GENOME in *_trimfilt5.fastq.gz
do 
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
    racon -m 8 -x -6 -g -8 -w 500 -t $THREADS ${GENOME} ${HEADER}_assemblymapped.sam ${HEADER}_assembly.fasta > ${HEADER}_racon1.fasta
done
conda deactivate

# Generate .sam file for 2nd round of Racon
conda activate minimap2_samtools

for GENOME in *_trimfilt5.fastq.gz
do 
    HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
    minimap2 -ax map-ont -t $THREADS ${HEADER}_racon1.fasta ${GENOME} > ${HEADER}_racon1mapped.sam
done
conda deactivate

# Run Racon for the 2nd time
conda activate racon
conda list

for GENOME in *_trimfilt5.fastq.gz
do 
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
    racon -m 8 -x -6 -g -8 -w 500 -t $THREADS ${GENOME} ${HEADER}_racon1mapped.sam ${HEADER}_racon1.fasta > ${HEADER}_racon2.fasta
done
conda deactivate

