###################################################################################
## This script polishes long-read assemblies with short reads using Polypolish
## (https://github.com/rrwick/Polypolish) DOI:10.1371/journal.pcbi.1009802. 
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh


# Run polypolish to polish long-read assemblies with short reads
conda activate polypolish_v0_5_0
conda list
mkdir 14.polypolish

for GENOME in *medakaconsensus.fasta
do 
	HEADER=$(echo ${GENOME} | sed 's/_flye00_medakaconsensus.fasta//');
	bwa index ${GENOME}
	bwa mem -t ${THREADS} -a ${GENOME} 12.illumina/${HEADER}_R1_illumina_cleantrim.fastq.gz > ${HEADER}_alignments_1.sam
	bwa mem -t ${THREADS} -a ${GENOME} 12.illumina/${HEADER}_R2_illumina_cleantrim.fastq.gz > ${HEADER}_alignments_2.sam
	polypolish_insert_filter.py --in1 ${HEADER}_alignments_1.sam --in2 ${HEADER}_alignments_2.sam --out1 ${HEADER}_filtered_1.sam --out2 ${HEADER}_filtered_2.sam
	polypolish ${GENOME} ${HEADER}_filtered_1.sam ${HEADER}_filtered_2.sam > ${HEADER}_polypolished.fasta
done

rm *.sam
rm *.amb 
rm *.ann 
rm *.bwt 
rm *.pac 
rm *.sa
mv *_medakaconsensus.fasta 14.polypolish/
echo "Polypolishing complete"
conda deactivate
