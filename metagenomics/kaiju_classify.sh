#!/bin/bash

#SBATCH --time=00:44:59
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=kaiju_classify
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --output=kaiju_class.log
#SBATCH --error=kaiju_class.err

module load StdEnv/2020 gcc/9.3.0 kaiju/1.7.4

kaiju -h

for R1 in /home/scratch/210610_Illumina/*R1_illumina_cleantrim.fastq.gz
do
        HEADER=$(echo ${R1} | sed 's/_R1_illumina_cleantrim.fastq.gz//');
	R2=$(echo ${R1} | sed 's/R1/R2/');
	kaiju -t nreuk_db_210728/nodes.dmp -f nreuk_db_210728/nr_euk/kaiju_db_nr_euk.fmi -i $R1 -j $R2 -z 32 -o ${HEADER}_kaiju.out
	kaiju-addTaxonNames -t nreuk_db_210728/nodes.dmp -n nreuk_db_210728/names.dmp -i ${HEADER}_kaiju.out -o ${HEADER}_kaiju_names.out
	kaiju2table -t nreuk_db_210728/nodes.dmp -n nreuk_db_210728/names.dmp -r phylum -m 1.0 -o ${HEADER}_kaiju_summary.tsv ${HEADER}_kaiju.out
done

echo "kaiju euk complete"
