###################################################################################
## This script was written to find the consensus sequence for long-read assembled
## samples that had reads depths greater than the target depth using Trycycler
## (https://github.com/rrwick/Trycycler) DOI:10.1186/s13059-021-02483-z
##
## Manual intervention is required when running Trycycler.
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8
TARGET_DEPTH=50      

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Use Trycycler to find the consensus sequences for samples that have read depths > $TARGET_DEPTH
conda activate trycycler_v0_5_3
conda list
trycycler --version

# Run trycycler
cd 05.seqtk/

for GENOME in *_trimfilt5.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
	trycycler cluster --threads $THREADS --assemblies ../10.trycycler/${HEADER}/*.fasta --reads ${GENOME} --out_dir ../10.trycycler/${HEADER}_trycycler 
done

cd ../
echo "Trycycler cluster complete"

# Manually inspect clusters for each samples, and remove those that are spurious
# rm -r cluster_00x

# Reconcile remaining clusters 
# trycycler reconcile --threads $THREADS --reads ../05.seqtk/${HEADER}_trimfilt5.fastq --cluster_dir ${HEADER}_trycycler/cluster_001

# Run a multiple sequence alignment
# trycycler msa --threads $THREADS --cluster_dir ${HEADER}_trycycler/cluster_001

# Partition reads into each contig
# trycycler partition --threads $THREADS --cluster_dirs ${HEADER}_trycycler/cluster_* --reads ../05.seqtk/${HEADER}_trimfilt5.fastq

# Generate consensus sequences for each contig
# trycycler consensus --threads $THREADS --cluster_dir ${HEADER}_trycycler/cluster_001

# Concatenate all contigs in a sample into one multi-fasta file
# cat ${HEADER}_trycycler/cluster_*/7_final_consensus.fasta > ${HEADER}_trycyclerconsensus.fasta

# Move to medaka folder
# mkdir ../13.medaka/
# mv ${HEADER}_trycycler/cluster_*/${HEADER}_trycyclerconsensus.fasta ../13.medaka/
