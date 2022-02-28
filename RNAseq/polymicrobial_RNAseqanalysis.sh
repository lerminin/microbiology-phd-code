###################################################################################
## This script processes multiple raw .fastq files from the sequencers all the way 
## through to gene counts and TPMs. Polymicrobial-specific commands come
## in bowtie2 and ensuring your concatenated cds input file has appropriate names.
##
## TPMs are generated with Salmon (https://salmon.readthedocs.io/en/latest/salmon.html)
###################################################################################

#!/bin/bash

# Must specify REFERENCE (name of transcriptome) for bowtie2 alignment
# Should have the same prefix as the concatenated fasta file
REFERENCE=SLandECWandplas
THREADS=8

# Remove phiX contamination from fastq files
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate bbtools
mkdir 01.bbduk.phiXcontam/
# If directory already exists, Unix will not overwrite it, but will overwrite files within if they have the same name
 
for GENOME in *_001.fastq.gz
 do
	HEADER=$(echo ${GENOME} | sed 's/_S.*_R1_.*.fastq.gz//');
	bbduk.sh in=${GENOME} out=${HEADER}_clean.fastq.gz outm=${HEADER}_phiXcontam.fastq.gz \
	ref=/opt/tools/bbmapResources/phix174_ill.ref.fa k=31 hdist=1 stats=${HEADER}_filteredphiXcontam.txt;
	mv *_phiXcontam.fastq.gz 01.bbduk.phiXcontam/;
	mv *_filteredphiXcontam.txt 01.bbduk.phiXcontam/
done

echo "phiX decontamination complete"
conda deactivate

# Run fastqc to check sequencing quality and to inform trimmomatic variables 
# Fastqc manual: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ 
conda activate fastqc
mkdir 02.fastqc.rawIlluminaReads/

for GENOME in *_clean.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_clean.fastq.gz//');
	fastqc ${GENOME};
	mv *fastqc.html 02.fastqc.rawIlluminaReads/;
	mv *fastqc.zip 02.fastqc.rawIlluminaReads/
done

echo "Fastqc analysis complete"
conda deactivate

# Run trimmomatic to remove adaptors and low quality base calls
# Trimmomatic manual: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf 
conda activate trimmomatic
mkdir 03.trimmomatic.trimmedIlluminaReads/

for GENOME in *_clean.fastq.gz
do
	HEADER=$(echo ${GENOME} | sed 's/_clean.fastq.gz//');
	trimmomatic SE -threads $THREADS -trimlog ${HEADER}_trimmomaticLogfile.log ${GENOME} \
	${HEADER}_cleanTrimmed.fastq \
	ILLUMINACLIP:/opt/tools/trimmomaticResources/illuminaAdaptors/TruSeq3-SE.fa:2:30:10 \
	CROP:73 HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36;
	mv *Logfile.log 03.trimmomatic.trimmedIlluminaReads/
done

echo "Trimmomatic analysis complete"
conda deactivate 

# Run fastqc to check sequencing quality after trimmomatic trimming
conda activate fastqc
mkdir 02.fastqc.trimmedIlluminaReads/

for GENOME in *_cleanTrimmed.fastq
do 
	HEADER=$(echo ${GENOME} | sed 's/_cleanTrimmed.fastq//');
	fastqc ${GENOME};
	mv *fastqc.html 02.fastqc.trimmedIlluminaReads/;
	mv *fastqc.zip 02.fastqc.trimmedIlluminaReads/
done

echo "Fastqc analysis complete"
conda deactivate

# Run bowtie2 to align reads to transcriptome file (cds)
# Bowtie2 manual: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# Make one .fa file with all coding sequences from both species = SLandECWandplas_cds_named.fa -> need REFERENCE name above
# cat *.fa > ${REFERENCE}_cds_named.fa
# When I tested this with bowtie2 v2.3.5, the input order of .fa had no impact on number of genes mapping
# By default, if a read aligns in more than one spot and has identical scores, bowtie2 will randomly pick one for the primary alignment
# BUT this can be changed to include secondary alignments 

conda activate bowtie2
mkdir 04.bowtie2.indexes/
mkdir 04.bowtie2.samfiles/

bowtie2-build ${REFERENCE}_cds_named.fa ${REFERENCE}

for GENOME in *_cleanTrimmed.fastq
do
	HEADER=$(echo ${GENOME} | sed 's/_cleanTrimmed.fastq//');
	bowtie2 -q -p $THREADS --phred33 -x ${REFERENCE} --very-sensitive --nofw \
	--un-gz /mnt/vault/personalFolders/nicole/polymicrobial_RNAseq_2020/un-seqs_${HEADER} \
	-U ${GENOME} -S ${HEADER}_aligned${REFERENCE}.sam;
	mv un-seqs_* 04.bowtie2.samfiles/
done

echo "Bowtie2 alignment complete"
conda deactivate

# Can check FLAGS in .sam files to see how reads are mapping 
# awk '{print $2}' filename.sam > FLAGS_filename.txt
# sort FLAGS_filename.txt | uniq -c | head -n 10

# Convert .sam to .bam files
# samtools manual: http://www.htslib.org/doc/samtools.html
conda activate samtools
mkdir 05.samtools

for GENOME in *.sam
do
	HEADER=$(echo ${GENOME} | sed 's/_aligned.*.sam//');
	samtools view -Sb -F 0x4 ${GENOME} > ${HEADER}.bam;
	samtools sort ${HEADER}.bam -o ${HEADER}_sorted.bam;
	samtools index ${HEADER}_sorted.bam;
	mv *_sorted.bam 05.samtools/;
	mv *.bai 05.samtools/
done

echo "Samtools conversion complete"
conda deactivate

# Calculate TPMs and counts with Salmon using bowtie2 alignments
# Salmon manual: https://salmon.readthedocs.io/en/latest/salmon.html
conda activate salmon

for GENOME in *.bam
do
	HEADER=$(echo ${GENOME} | sed 's/.bam//');
	salmon quant -t ${REFERENCE}_cds_named.fa -l SR -a ${GENOME} -o ${HEADER}_alignmode
done

echo "Salmon analysis complete"
conda deactivate 

