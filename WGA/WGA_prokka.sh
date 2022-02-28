###################################################################################
## This script takes the completed and assembled genome and does some initial 
## annotation and investigation:
##
## Annotates using Prokka (https://github.com/tseemann/prokka) 
## Plasmid detection by MOB-suite (https://github.com/phac-nml/mob-suite)
## Identification of the following through Abricate (https://github.com/tseemann/abricate) 
##  - Plasmid replicon detection by PlasmidFinder (https://cge.cbs.dtu.dk/services/PlasmidFinder/) 
##  - Antimicrobial resistance genes in CARD (https://card.mcmaster.ca/)
##  - Antimicrobial resistance genes in ResFinder (https://cge.cbs.dtu.dk/services/ResFinder/) 
###################################################################################

#!/bin/bash

# Set number of cores, target depth for subsampling, and genome size parameters
THREADS=8

# Recognize conda
source /opt/miniconda3/etc/profile.d/conda.sh

# Annotate with prokka
conda activate prokka
conda list
mkdir 15.prokka/

for GENOME in *.fasta
do
	HEADER=$(echo ${GENOME} | sed 's/_polca.fasta//');
	prokka --cpus $THREADS --outdir ${HEADER}_polca ${GENOME};
	cp 15.prokka/${HEADER}_polca/PROKKA*.gff 15.prokka/${HEADER}_polca/${HEADER}_polca_prokkaannot.gff;
	mv 15.prokka/${HEADER}_polca/${HEADER}_polca_prokkaannot.gff . 
done

echo "Prokka annotation complete"
conda deactivate

# Identify plasmids with mobsuite
conda activate mobsuite_v3_0_0
conda list
mkdir 16.mobsuite

for GENOME in *.fasta
do
	HEADER=$(echo ${GENOME} | sed 's/_polca.fasta//');
	mob_typer --multi --infile ${GENOME} -o 16.mobsuite/${HEADER}_mobsuite.txt -n $THREADS
done

echo "Mobsuite classification complete"
conda deactivate

# Identify plasmid replicons and ARGs with abricate (plasmidfinder, CARD, ResFinder)
conda activate abricate
conda list
abricate --list
mkdir 17.abricate

for GENOME in *.fasta
do
	HEADER=$(echo ${GENOME} | sed 's/_polca.fasta//');
	abricate --threads $THREADS --db plasmidfinder ${GENOME} > ${HEADER}_polca_plasmidfinder.txt
    abricate --threads $THREADS --db resfinder ${GENOME} > ${HEADER}_polca_resfinder.txt
    abricate --threads $THREADS --db card ${GENOME} > ${HEADER}_polca_card.txt
done

echo "Abricate annotation complete"
conda deactivate

# Predict viral contigs with DeepVirFinder
conda activate deepvirfinder
conda list
mkdir 18.deepvirfinder

for GENOME in *.fasta
do
	HEADER=$(echo ${GENOME} | sed 's/_polca.fasta//');
	mkdir 18.deepvirfinder/${HEADER}
	python dvf.py -i ${GENOME} -o 18.deepvirfinder/${HEADER}/ -l 300 -c $THREADS
done

echo "DeepVirFinder prediction complete"
conda deactivate
