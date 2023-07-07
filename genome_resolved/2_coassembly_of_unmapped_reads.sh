#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=420gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/unmapped_reads_to_bins/unmapped_reads

megahit -1 /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/unmapped_reads_to_bins/unmapped_reads/all_unmapped_reads_R1.fastq -2 /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/unmapped_reads_to_bins/unmapped_reads/all_unmapped_reads_R2.fastq --k-min 31 --k-max 121 --k-step 10 --mem-flag 1 -m 429496729600 -t 20 -o unmapped_reads_coassembly_megahit_out
