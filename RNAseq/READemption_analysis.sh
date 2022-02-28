#############################################
# READemption analysis for polymicrobial RNAseq
# Assumes that: 
# 1. READemption 'create' command has been run 
# 2. All input files are in their folders
# 3. The project path folder has been specified
# 
# https://reademption.readthedocs.io/en/latest/
#############################################
#!/usr/bin/bash

# Set number of threads
THREADS=8

# Activate environment
source /opt/tools/virtualenv/READemption/bin/activate
reademption -h

# align
reademption align -p $THREADS -g --fastq --segemehl_evalue 0.001 --project_path READemption_analysis/ --crossalign_cleaning "Salmonella:NC_017718.1,NC_017719.1,NC_017720.1,NC_016810.1;Ecoli:NC_007779.1,AJ310442.1"

# coverage 
reademption coverage -p $THREADS --project_path READemption_analysis/

# gene_quanti 
reademption gene_quanti -p $THREADS --project_path READemption_analysis/

# deseq
reademption deseq -l WT_WT_7A_cleanTrimmed.fastq,WT_WT_7B_cleanTrimmed.fastq,cib_WT_7A_cleanTrimmed.fastq,cib_WT_7B_cleanTrimmed.fastq,WT_cirA_7A_cleanTrimmed.fastq,WT_cirA_7B_cleanTrimmed.fastq -c WT_WT,WT_WT,cib_WT,cib_WT,WT_cirA,WT_cirA --project_path READemption_analysis/

# viz-align 
reademption viz_align --project_path READemption_analysis/

# viz-gene-quanti
reademption viz_gene_quanti --project_path READemption_analysis/

# viz-deseq
reademption viz_deseq --project_path READemption_analysis
