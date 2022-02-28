###################################################################################
## This script polishes long-read assemblies with short reads using POLCA 
## (https://github.com/alekseyzimin/masurca) DOI:10.1371/journal.pcbi.1007981
##
## To run, need to specify specific short read names in polca.sh command
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Run polca to polish long-read assemblies with short reads - cannot loop, need to change polca command with specific short read names
conda activate polca_v4_0_7
conda list
mkdir 15.polca

for GENOME in *polypolished.fasta
do 
	HEADER=$(echo ${GENOME} | sed 's/_polypolished.fasta//');
	polca.sh -a ${GENOME} -r 'barXX_R1_illumina_cleantrim.fastq.gz barXX_R2_illumina_cleantrim.fastq.gz' -t ${THREADS} ## add short read names
	mv *.report 15.polca/
	mv ${GENOME}.PolcaCorrected.fa ${HEADER}_polca.fasta
	mkdir 15.polca/${HEADER}
	mv *fasta.* 15.polca/${HEADER}
done

rm *.err
mv *_polypolished.fasta 15.polca
echo "POLCA complete"
conda deactivate
