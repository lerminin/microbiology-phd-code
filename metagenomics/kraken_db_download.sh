#!/bin/bash

#SBATCH --time=167:59:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=kraken_db
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --output=kraken_db.log
#SBATCH --error=kraken_db.err

module load StdEnv/2020 gcc/9.3.0 kraken2/2.1.1

# Download taxonomy - took < 1 h
#kraken2-build --download-taxonomy --db nt_210709

# Download taxa - took 10 h and <25 G memory for all
#kraken2-build --download-library bacteria --use-ftp --db bac_vir_arch_proto_plant_210707
#kraken2-build --download-library archaea --use-ftp --db bac_vir_arch_proto_plant_210707
#kraken2-build --download-library viral --use-ftp --db bac_vir_arch_proto_plant_210707
#kraken2-build --download-library protozoa --use-ftp --db bac_vir_arch_proto_plant_210707
#kraken2-build --download-library plant --use-ftp --db bac_vir_arch_proto_plant_210707
#kraken2-build --download-library plasmid --use-ftp --db bac_vir_arch_proto_plant_210707
#kraken2-build --download-library fungi --use-ftp --db bac_vir_arch_proto_plant_210707
# Took 42 h for nt database
#kraken2-build --download-library nt --use-ftp --db nt_210709

# Build database - took 32 threads, 14 h and <128 G memory for bac_vir_arch_proto_plant_210707
kraken2-build --build --db nt_210709 --threads 32
