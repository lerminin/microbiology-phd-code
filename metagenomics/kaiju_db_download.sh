#!/bin/bash

#SBATCH --time=17:59:00
#SBATCH --mail-user=email@email.email
#SBATCH --mail-type=ALL
#SBATCH --job-name=kaiju_db
#SBATCH --cpus-per-task=16
#SBATCH --mem=450G
#SBATCH --output=kaiju_db.log
#SBATCH --error=kaiju_db.err

module load StdEnv/2020 gcc/9.3.0 kaiju/1.7.4

kaiju -h 

# Took about 9 hours and 400 G memory with 16 threads for nr
# Took about 8 hours and 450 G memory with 16 threads for nr_euk
mkdir nreuk_db_210729 
cd nr_db_210729
kaiju-makedb -s nr_euk -t 16
