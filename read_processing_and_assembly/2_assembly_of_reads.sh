#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=420gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads

for element in $(<$1)
do 
megahit -1 "$element"_1_trimmed.fastq -2 "$element"_2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o "$element"_megahit_out
done
