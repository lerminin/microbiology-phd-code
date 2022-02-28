# microbiology-phd-code
Code for Microbiology PhD Thesis - bioinformatic whole genome assembly, RNAseq data processing and analysis, 16S analysis, figures in R

Scripts are divided into the following sections:

1. [Whole genome assembly (WGA)](#whole-genome-assembly-wga-folder)
2. [RNAseq](#rnaseq-folder)
3. [Metagenomics](#metagenomics-folder)
4. [Conda environments](#conda-environments-conda_envs-folder) 

## Whole genome assembly (WGA) folder

These scripts were used to assemble whole genomes for individual bacterial isolates which were sequenced with both Illumina paired-end short reads and Oxford Nanopore Technologies long reads. The scripts are designed to run in the order specified below and are constructed to loop through multiple samples with minimal intervention required. Quality checks are included throughout. 

### guppy_basecalling.sh

This script was run on a computing cluster to basecall Oxford Nanopore long reads using the SUP model for MinION data using GPUs. 

### WGA_readqc_assemblies.sh

This script processes multiple raw .fastq files from the sequencers all the way through to assembly and annotation, and includes the following programs:

- Barcodes are removed with [Porechop](https://github.com/rrwick/Porechop)
- Reads are trimmed with [Filtlong](https://github.com/rrwick/Filtlong) 
- Read quality metrics are assessed with [Nanopack](https://github.com/wdecoster/nanopack) ([DOI:10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149)) 
- Read subsets generated with [seqtk](https://github.com/lh3/seqtk)
- Assemblies by [Flye](https://github.com/fenderglass/Flye) ([DOI:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)) 
- Assemblies by [Raven](https://github.com/lbcb-sci/raven) ([DOI:10.1038/s43588-021-00073-4](http://dx.doi.org/10.1038/s43588-021-00073-4))
- Assemblies by [Canu](https://github.com/marbl/canu) ([DOI:10.1101/gr.215087.116](https://doi.org/10.1101/gr.215087.116))
- Assemblies by [Redbean/wtdbg2](https://github.com/ruanjue/wtdbg2) ([DOI:10.1038/s41592-019-0669-3](https://www.nature.com/articles/s41592-019-0669-3)) 
- Assemblies by [Miniasm/minipolish](https://github.com/lh3/miniasm) ([DOI:10.1093/bioinformatics/btw152](https://doi.org/10.1093/bioinformatics/btw152))  

Prior to running this script, all fastq files must be concatenated into one .fastq.gz per barcode and be in the working directory. This script is constructed so that if a sample has fewer reads than the required depth for Trycycler (at least 50X), it will not be subsampled and will be assembled by Flye.

### WGA_trycycler.sh 

This script finds the consensus sequence between different long-read assemblies for samples that had reads depths greater than the target depth using [Trycycler](https://github.com/rrwick/Trycycler) ([DOI:10.1186/s13059-021-02483-z](https://doi.org/10.1186/s13059-021-02483-z)). Assemblies are generated by the WGA_readqc_assemblies.sh script. Manual intervention and decision making is required when running Trycycler.

### WGA_medaka.sh 

This script polishes Oxford Nanopore long-read assemblies from Trycycler with long reads using [Medaka](https://github.com/nanoporetech/medaka) and CPUs.

### WGA_illumina.sh

The script processes raw multiple paired-end Illumina reads, and includes the following programs:

- Barcodes are removed with [BBDuk](https://sourceforge.net/projects/bbmap/)
- Reads are trimmed with [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) ([DOI:10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)) 
- Read quality metrics are assessed with [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### WGA_unicycler.sh

This script generates hybrid assemblies using Illumina paired-end short reads and Oxford Nanopore long reads with [Unicycler](https://github.com/rrwick/Unicycler) ([DOI:10.1371/journal.pcbi.1005595](https://doi.org/10.1371/journal.pcbi.1005595)). 

### WGA_polypolish.sh

This script polishes long-read assemblies with short reads using [Polypolish](https://github.com/rrwick/Polypolish) ([DOI:10.1371/journal.pcbi.1009802](https://doi.org/10.1371/journal.pcbi.1009802). This script runs on assemblies that have already been polished by Medaka using WGA_medaka.sh. 

### WGA_polca.sh

This script polishes long-read assemblies with short reads using [POLCA](https://github.com/alekseyzimin/masurca) ([DOI:10.1371/journal.pcbi.1007981](https://doi.org/10.1371/journal.pcbi.1007981). This script runs on assemblies that have already been polished by PolyPolish using WGA_polypolish.sh. 

### WGA_pilon.sh

After Medaka long-read polishing, an alternative to Polypolish & POLCA for polishing assemblies with Illumina reads is [Pilon](https://github.com/broadinstitute/pilon) ([DOI:10.1371/journal.pone.0112963](https://doi.org/10.1371/journal.pone.0112963)). The script makes use of the `insertsizeI.py` and `insertsizeX.py` scripts for determining the minimum and maximum insert size between paired end reads for the bowtie2 call. 

### WGA_racon.sh

After Medaka long-read polishing, an additional long-read polishing step can be done with [Racon](https://github.com/lbcb-sci/racon) ([DOI:10.1101/gr.214270.116](https://dx.doi.org/10.1101%2Fgr.214270.116)). 

### WGA_prokka.sh

This script takes the completed and assembled genome and does some initial annotation and investigation:

- Annotates using [Prokka](https://github.com/tseemann/prokka) ([DOI:10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153))
- Plasmid detection by [MOB-suite](https://github.com/phac-nml/mob-suite) ([DOI:10.1099/mgen.0.000206](https://doi.org/10.1099/mgen.0.000206))
- Identification of the following through [Abricate](https://github.com/tseemann/abricate)
  - Plasmid replicon detection by [PlasmidFinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/) ([DOI:10.1128/AAC.02412-14](https://doi.org/10.1128/aac.02412-14))
  - Antimicrobial resistance genes in [CARD](https://card.mcmaster.ca/) ([DOI:10.1093/nar/gkw1004](https://doi.org/10.1093/nar/gkw1004))
  - Antimicrobial resistance genes in [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) ([DOI:10.1093/jac/dkaa345](https://doi.org/10.1093/jac/dkaa345))
- Prediction of viral contigs using [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) ([DOI:10.1007/s40484-019-0187-4](https://doi.org/10.1007/s40484-019-0187-4))

## RNAseq folder 

This folder contains scripts for processing and quantifying transcriptomes of a mixed-species bacterial culture generated from Illumina unpaired read sequencing data.  

### READemption_analysis.sh

This script runs the analysis pipeline [READemption](https://reademption.readthedocs.io/en/latest/) ([DOI:10.1093/bioinformatics/btu533](https://doi.org/10.1093/bioinformatics/btu533)) for quantifying RNAseq transcripts and assumes that the `create` command has already been run and that the input files are in the folder structure specified by the READemption documentation.

In addition, this script shows how to make use of the cross-align clean option to eliminate reads that map to multiple replicons during the `align` step. 

### polymicrobial_RNAseqanalysis.sh

This script runs the analysis pipeline for quantifying RNAseq transcripts using the following programs:

- Reads aligned to reference with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) ([DOI:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923))
- .sam files generated with [Samtools](http://www.htslib.org/doc/) ([DOI:10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352))
- TPMs and counts quantified with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) ([DOI:10.1038/nmeth.4197](https://dx.doi.org/10.1038%2Fnmeth.4197))

### DESeq2_analysis.R

This script runs the the analysis for differentially expressed genes using [DESeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) ([DOI:10.18129/B9.bioc.DESeq2](https://doi.org/doi:10.18129/B9.bioc.DESeq2)) and produces diagnostic and results plots. 

## Metagenomics folder

This folder contains scripts to assemble and investigate metagenomic DNA samples on a high performance computing cluster. 

### megahit.sh

This script assembles metagenomic Illumina paired end reads using [MEGAHIT](https://github.com/voutcn/megahit) ([DOI:10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)) with preset `meta-large`. 

### pspades.sh

This script assembles metagenomic Illumina paired end reads using [SPAdes](https://github.com/ablab/spades) ([DOI:10.1093/bioinformatics/btw493](https://doi.org/10.1093/bioinformatics/btw493)) with arguments `--meta` and `--plasmid`. 

### idba.sh

This script assembles metagenomic Illumina paired end reads using [IDBA-UD](https://github.com/loneknightpy/idba) ([DOI:10.1093/bioinformatics/bts174](https://doi.org/10.1093/bioinformatics/bts174)). 

### kaiju_db_download.sh and kaiju_classify.sh

This script downloads the database and classifies reads into taxonomic groups using [Kaiju](https://github.com/bioinformatics-centre/kaiju) ([DOI:10.1038/ncomms11257](https://www.nature.com/articles/ncomms11257)). 

### kraken_classify.sh and kraken_db_download.sh

This script downloads the database and classifies reads into taxonomic groups using [Kraken2](https://github.com/DerrickWood/kraken2) ([DOI:10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)). 

### metaxa2.sh

This script identifies the proportion of reads with small subunit and large subunit rRNA sequences using [METAXA2](https://microbiology.se/software/metaxa2/) ([DOI:10.1111/1755-0998.12399](https://doi.org/10.1111/1755-0998.12399)). 

## Conda environments (conda_envs) folder

Many of the shell scripts included here are dependent on conda environments run on a Linux Ubuntu x86_64 machine. The packages and versions for each environment used in the respective folders can be installed using the .yml files. 

