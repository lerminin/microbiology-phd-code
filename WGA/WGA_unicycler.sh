###################################################################################
## This script was written to generate hybrid assemblies using Unicycler 
## (https://github.com/rrwick/Unicycler) DOI:10.1371/journal.pcbi.1005595
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Use unicycler to do hybrid assembly
conda activate unicycler_v5_0
conda list
mkdir 11.unicycler
cd 11.unicycler 

for GENOME in *.fastq.gz
do 
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
	unicycler -1 ../12.illumina/${HEADER}_R1_illumina_cleantrim.fastq.gz -2 ../12.illumina/${HEADER}_R2_illumina_cleantrim.fastq.gz -s ../12.illumina/${HEADER}_U12_illumina_cleantrim.fastq.gz \
	-l ${HEADER}_trimfilt5.fastq.gz -o ${HEADER} -t $THREADS
done

cd ../
echo "Unicycler complete"
conda deactivate
