###################################################################################
## This script polishes Medaka-polished long read assemblies with paired-end Illumina
## reads using Pilon (https://github.com/broadinstitute/pilon). 
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=12

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Run pilon to polish draft assemblies with Illumina reads
conda activate pilonComplete
conda list 
mkdir 14.pilon

# First, build initial bowtie2 index
for GENOME in *.fasta
do 
	HEADER=$(echo ${GENOME} | sed 's/_medakaconsensus.fasta//');
	bowtie2-build ${GENOME} 14.pilon/${HEADER}_medakaconsensus;
	bowtie2 -1 12.illumina/${HEADER}_R1_illumina_cleantrim.fastq.gz -2 12.illumina/${HEADER}_R2_illumina_cleantrim.fastq.gz \
	-x 14.pilon/${HEADER}_medakaconsensus --sensitive --threads $THREADS -X 1000 -S insertsize.sam
done

# Pass python variables to bash:
INSERTI=$(python3 insertsizeI.py)
echo $INSERTI
INSERTX=$(python3 insertsizeX.py)
echo $INSERTX

# Generate first bam files and do initial Pilon run
for GENOME in *.fasta
do 
	HEADER=$(echo ${GENOME} | sed 's/_medakaconsensus.fasta//');
	bowtie2 -1 12.illumina/${HEADER}_R1_illumina_cleantrim.fastq.gz -2 12.illumina/${HEADER}_R2_illumina_cleantrim.fastq.gz \
	-x 14.pilon/${HEADER}_medakaconsensus --very-sensitive-local --threads $THREADS -I $INSERTI -X $INSERTX | samtools sort > 14.pilon/${HEADER}_illuminamapped.bam;
	bowtie2 -U 12.illumina/${HEADER}_U12_illumina_cleantrim.fastq.gz -x 14.pilon/${HEADER}_medakaconsensus --very-sensitive-local --threads $THREADS \
	| samtools sort > 14.pilon/${HEADER}_illuminamappedU.bam;
	samtools index 14.pilon/${HEADER}_illuminamapped.bam;
	samtools index 14.pilon/${HEADER}_illuminamappedU.bam;
	pilon --genome ${GENOME} --frags 14.pilon/${HEADER}_illuminamapped.bam --unpaired 14.pilon/${HEADER}_illuminamappedU.bam --output ${HEADER}_round1 --changes --threads $THREADS;
	head -n 20 ${HEADER}_round1.changes;
	sed -i 's/_pilon//' ${HEADER}_round1.fasta 
	mv ${GENOME} 14.pilon/
done

# Run Pilon for another few rounds
for GENOME in *.fasta
do
        for i in {1..6}; do # can adjust numbers to do more rounds, i.e. {1..10}
                HEADER=$(echo ${GENOME} | sed 's/_round1.fasta//');
                x=$((i+1));
                bowtie2-build ${HEADER}_round$i.fasta 14.pilon/${HEADER}_round$i;
                bowtie2 -1 12.illumina/${HEADER}_R1_illumina_cleantrim.fastq.gz -2 12.illumina/${HEADER}_R2_illumina_cleantrim.fastq.gz \
                -x 14.pilon/${HEADER}_round$i --very-sensitive-local --threads $THREADS -I $INSERTI -X $INSERTX | samtools sort > 14.pilon/${HEADER}_illuminamapped.bam;
                bowtie2 -U 12.illumina/${HEADER}_U12_illumina_cleantrim.fastq.gz -x 14.pilon/${HEADER}_round$i --very-sensitive-local --threads $THREADS \
                | samtools sort > 14.pilon/${HEADER}_illuminamappedU.bam;
                samtools index 14.pilon/${HEADER}_illuminamapped.bam;
                samtools index 14.pilon/${HEADER}_illuminamappedU.bam;
                pilon --genome ${HEADER}_round$i.fasta --frags 14.pilon/${HEADER}_illuminamapped.bam --unpaired 14.pilon/${HEADER}_illuminamappedU.bam --output ${HEADER}_round$x --changes --threads $THREADS;
                head -n 20 ${HEADER}_round$x.changes;
                sed -i 's/_pilon//' ${HEADER}_round$x.fasta
        done
done

# Examine .changes files and see which round*.fasta file to keep
#cp *round*.fasta *pilonconsensus.fasta

mv *round* 14.pilon/
rm 14.pilon/*.bt2
echo "Pilon polishing complete"
conda deactivate
