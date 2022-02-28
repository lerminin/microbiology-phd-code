# microbiology-phd-code
Code for Microbiology PhD Thesis - bioinformatic whole genome assembly, 16S analysis, figures in R


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

### WGA_prokka.sh

This script takes the completed and assembled genome and does some initial annotation and investigation:

- Annotates using [Prokka](https://github.com/tseemann/prokka) ([DOI:10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153))
- Plasmid detection by [MOB-suite](https://github.com/phac-nml/mob-suite) ([DOI:10.1099/mgen.0.000206](https://doi.org/10.1099/mgen.0.000206))
- Identification of the following through [Abricate](https://github.com/tseemann/abricate)
  - Plasmid replicon detection by [PlasmidFinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/) ([DOI:10.1128/AAC.02412-14](https://doi.org/10.1128/aac.02412-14))
  - Antimicrobial resistance genes in [CARD](https://card.mcmaster.ca/) ([DOI:10.1093/nar/gkw1004](https://doi.org/10.1093/nar/gkw1004))
  - Antimicrobial resistance genes in [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) ([DOI:10.1093/jac/dkaa345](https://doi.org/10.1093/jac/dkaa345))

## Conda environments (conda_envs) folder

The scripts included here are dependent on conda environments run on a Linux Ubuntu x86_64 machine. The packages and versions for each environment used in the respective folders can be installed using the .yml files. 

