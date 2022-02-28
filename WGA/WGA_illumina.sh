###################################################################################
## This script was written to process multiple paired-end Illumina .fastq files 
## and includes quality checks, barcode trimming, and filtering.
## 
## Rename illumina reads to HEADER_R*_illumina.fastq.gz and put in illumina/ folder
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Remove phiX contamination from Illumina reads
# Rename illumina reads to HEADER_R*_illumina.fastq.gz and put in illumina/ folder
conda activate bbtools
conda list 
mkdir 12.illumina/

for GENOME in *_trimfilt5
do
	cd 12.illumina/
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5//');
	bbduk.sh in=${HEADER}_R1_illumina.fastq.gz in2=${HEADER}_R2_illumina.fastq.gz out=${HEADER}_R1_illumina_clean.fastq.gz out2=${HEADER}_R2_illumina_clean.fastq.gz \
	outm=${HEADER}_R1_illumina_phiXcontam.fastq.gz outm2=${HEADER}_R2_illumina_phiXcontam.fastq.gz ref=/opt/tools/bbmapResources/phix174_ill.ref.fa k=31 hdist=1 \
	stats=${HEADER}_bbduk_filteredphiXcontam.txt
	cd ../
done

echo "Bbduk decontamination complete"
conda deactivate

# Use fastqc to check quality scores
conda activate fastqc
conda list 

for GENOME in *_trimfilt5
do
	cd 12.illumina/
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5//');
	fastqc -t ${THREADS} ${HEADER}_R1_illumina_clean.fastq.gz;
	fastqc -t ${THREADS} ${HEADER}_R2_illumina_clean.fastq.gz;
	rm *.zip
	cd ../
done

echo "Fastqc check complete"
conda deactivate

#Use trimmomatic to trim Illumina reads:
conda activate trimmomatic

for GENOME in *_trimfilt5
do
	cd 12.illumina/
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5//');
	trimmomatic PE -threads $THREADS -trimlog ${HEADER}_trimmomatic.log ${HEADER}_R1_illumina_clean.fastq.gz ${HEADER}_R2_illumina_clean.fastq.gz \
	${HEADER}_R1_illumina_cleantrim.fastq.gz ${HEADER}_U1_illumina_cleantrim.fastq.gz ${HEADER}_R2_illumina_cleantrim.fastq.gz ${HEADER}_U2_illumina_cleantrim.fastq.gz \
	ILLUMINACLIP:/opt/tools/trimmomaticResources/illuminaAdaptors/TruSeq3-SE.fa:2:30:10:8:TRUE CROP:270 HEADCROP:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
	cat ${HEADER}_U*_illumina_cleantrim.fastq.gz > ${HEADER}_U12_illumina_cleantrim.fastq.gz
	cd ../
done 

# Use fastqc to check quality scores
conda activate fastqc

for GENOME in *_trimfilt5
do
	cd 12.illumina/
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5//');
	fastqc -t ${THREADS} ${HEADER}_R1_illumina_cleantrim.fastq.gz;
	fastqc -t ${THREADS} ${HEADER}_R2_illumina_cleantrim.fastq.gz;
	rm *.zip
	cd ../
done

echo "Fastqc check complete"
conda deactivate

