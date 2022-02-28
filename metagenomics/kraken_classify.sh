#!/bin/bash

#SBATCH --time=00:59:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=kraken_classify
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=kraken_class.log
#SBATCH --error=kraken_class.err

module load StdEnv/2020 gcc/9.3.0 kraken2/2.1.1

kraken2 -h

# Run kraken2
kraken2 --paired --classified-out sample_class#.fastq ../210610_Illumina/sample_R1_illumina_cleantrim.fastq.gz ../210610_Illumina/sample_R2_illumina_cleantrim.fastq.gz --threads 32 \
	--gzip-compressed --db bac_vir_arch_proto_plant_210707 --output sample_profiled_metagenome.txt --report sample_kraken_report.txt --use-names

echo "kraken2 complete"
