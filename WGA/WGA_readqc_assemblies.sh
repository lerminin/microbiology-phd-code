###################################################################################
## This script was written to process multiple raw .fastq files from the sequencers 
## all the way through to assembly and annotation. Prior to running this script, 
## all fastq files must be concatenated into one .fastq.gz per barcode and be in  
## the working directory.
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=14
TARGET_DEPTH=50                 # 50x depth of coverage.
GENOME_SIZE=5200000
GENOME_SIZE_ASSEMBLY="5.2m"

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Remove barcodes
conda activate porechop
conda list
mkdir 01.porechop/
# If directory already exists, Unix will not overwrite it, but will overwrite files within if they have the same name
 
for GENOME in *.fastq.gz
 do
	HEADER=$(echo ${GENOME} | sed 's/.fastq.gz//');
	porechop -i ${GENOME} -o ${HEADER}_trim.fastq.gz -t $THREADS;
	mv ${GENOME} 01.porechop/
done

echo "Porechop complete"
conda deactivate

# Run nanoplot to check sequencing quality and to inform filtlong variables 
# Nanoplot manual: https://github.com/wdecoster/NanoPlot
conda activate nanopack
conda list
mkdir 02.nanoplot/

for GENOME in *_trim.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_trim.fastq.gz//');
	NanoPlot -t $THREADS --fastq ${GENOME} -o ${HEADER} --loglength;
	mv ${HEADER} 02.nanoplot/
done

echo "Nanoplot analysis complete"
conda deactivate

# Run filtlong to remove short reads and low quality base calls
# Filtlong manual: https://github.com/rrwick/Filtlong
# Filter settings:
#	filt5: --min_length 1000 --min_mean_q 90
conda activate filtlong
conda list
mkdir 03.filtlong/

for GENOME in *_trim.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_trim.fastq.gz//');
	filtlong --min_length 1000 --min_mean_q 90 ${GENOME} | gzip > ${HEADER}_trimfilt5.fastq.gz;
	mv ${GENOME} 03.filtlong/
done

echo "Filtlong analysis complete"
conda deactivate 

# Run nanoplot to check sequencing quality after trimmomatic trimming
conda activate nanopack
conda list
mkdir 04.nanoplot/

for GENOME in *_trimfilt5.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
	NanoPlot -t $THREADS --fastq ${GENOME} -o ${HEADER}_trimfilt5 --loglength;
	mv ${HEADER} 04.nanoplot/
done

echo "Nanoplot analysis complete"
conda deactivate

# Run seqtk to subsample reads. Set to 50X coverage and estimated genome size
conda activate seqtk
conda list 

mkdir 05.seqtk
mkdir 06.flye                     
mkdir 07.raven
mkdir 08.canu
mkdir 09.miniasm_minipolish
mkdir 10.trycycler
# If directory already exists, Unix will not overwrite it, but will overwrite files within if they have the same name


for GENOME in *_trimfilt5.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
	MEAN_READ_LENGTH=$(seqtk comp ${GENOME} | awk '{count++; bases += $2} END{print bases/count}');
	NUM_READS_FOR_TARGET_DEPTH=$(echo $TARGET_DEPTH"*"$GENOME_SIZE"/"$MEAN_READ_LENGTH | bc);
	NUM_READS_TOTAL=$(echo $(cat ${GENOME}|wc -l)/4 | bc);
	echo "$MEAN_READ_LENGTH";             # should match Nanoplot results
	echo "$NUM_READS_FOR_TARGET_DEPTH";
	echo "$NUM_READS_TOTAL";

	if (("$NUM_READS_FOR_TARGET_DEPTH" < "$NUM_READS_TOTAL"))
	then	
		for i in {0..11}; do
			x=""
			if (($i < 10))
				then x="0"
			fi

			seqtk sample -s $i ${GENOME} "$NUM_READS_FOR_TARGET_DEPTH" | paste - - - - | shuf | tr '\t' '\n' > ${HEADER}_sample"$x$i".fastq;
				
			# Move subsamples to their respective folders
			if (( $i == 0 || $i == 1 || $i == 2))
				then mv ${HEADER}_sample"$x$i".fastq 06.flye/
			elif (( $i == 3  || $i == 4  || $i == 5 ))
				then mv ${HEADER}_sample"$x$i".fastq 07.raven/;
			elif (( $i == 6  || $i == 7 || $i == 8 ))
				then mv ${HEADER}_sample"$x$i".fastq 08.canu/;
			else 
				mv ${HEADER}_sample"$x$i".fastq 09.miniasm_minipolish;
			fi
		done
	else
		mv ${GENOME} 06.flye
	fi

	mv ${GENOME} 05.seqtk/
done

# Use flye to assemble first set of subsamples
conda activate flye_v2_9
conda list
cd 06.flye/

for GENOME in *.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/.fastq//');
	flye --nano-hq ${GENOME} --out-dir ${HEADER} --threads $THREADS --genome-size ${GENOME_SIZE_ASSEMBLY};
	mv ${HEADER}/assembly.fasta ${HEADER}/${HEADER}_flye.fasta;
	mv ${HEADER}/${HEADER}_flye.fasta ../10.trycycler/
done

cd ../
echo "Flye assembly complete"
conda deactivate

# Use raven to assemble the second set of subsamples
conda activate raven_v1_7_0
conda list
cd 07.raven/

for GENOME in *.fastq
do
	HEADER=$(echo ${GENOME} | sed 's/.fastq//');
	raven -t $THREADS ${GENOME} > ${HEADER}_raven.fasta;
	mv ${HEADER}_raven.fasta ../10.trycycler/
done

rm raven.cereal
cd ../
echo "Raven assembly complete"
conda deactivate

# Use canu to assemble the third set of subsamples
# Note that canu scales itself accordingly to the number of cores available and there's no option to set threads, so use with caution
conda activate canu_v2_1_1
conda list
cd 08.canu/

for GENOME in *.fastq
do
	HEADER=$(echo ${GENOME} | sed 's/.fastq//');
	canu -p ${HEADER} -d ${HEADER} genomeSize=${GENOME_SIZE_ASSEMBLY} -nanopore ${GENOME};
	mv ${HEADER}/${HEADER}.contigs.fasta ../10.trycycler/
done

cd ../
echo "Canu assembly complete"
conda deactivate

## Alternative to canu, use redbean/wtdbg2 to assemble the third set of subsamples
#conda activate wtdbg_v2_5
#conda list
#cd 08.canu/

#for GENOME in *.fastq
#do
#	HEADER=$(echo ${GENOME} | sed 's/.fastq//');
#	wtdbg2 -x ont -t $THREADS -g ${GENOME_SIZE_ASSEMBLY} -i ${GENOME} -fo ${HEADER}_wtdbg;
#	wtpoa-cns -t $THREADS -i ${HEADER}_wtdbg.ctg.lay.gz -fo ${HEADER}_wtdbg.ctg.fasta;
#	mv ${HEADER}_wtdbg.ctg.fasta ../10.trycycler/
#done

#cd ../
#echo "Redbean/wtdbg2 assembly complete"
#conda deactivate

# Use miniasm/minipolish to assemble the fourth set of subsamples
conda activate mini_assembly
conda list
cd 09.miniasm_minipolish/

for GENOME in *.fastq
do
	HEADER=$(echo ${GENOME} | sed 's/.fastq//');
	minimap2 -t $THREADS -x ava-ont ${GENOME} ${GENOME} > ${HEADER}_overlaps.paf;
	miniasm -f ${GENOME} ${HEADER}_overlaps.paf > ${HEADER}_assembly.gfa;
	minipolish -t $THREADS ${GENOME} ${HEADER}_assembly.gfa > ${HEADER}_polishedassembly.gfa;                
	any2fasta ${HEADER}_polishedassembly.gfa > ${HEADER}_mini.fasta;
	mv ${HEADER}_mini.fasta ../10.trycycler/
done

cd ../
echo "Miniasm/minipolish assembly complete"
conda deactivate

# Sort samples into folders for trycycler prep
cd 05.seqtk
for GENOME in *_trimfilt5.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_trimfilt5.fastq.gz//');
	mkdir ../10.trycycler/${HEADER} 
done

cd ../10.trycycler

for GENOME in *.fasta
do
        HEADER=$(echo ${GENOME} | sed 's/_sample.*.fasta//');
        mv ${GENOME} ${HEADER}
done
cd ../


