#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies

for element in $(<$1)
do
bbmap.sh -Xmx48G threads=20 semiperfectmode=t overwrite=t ref=/home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/single_assembly_bins/all_255_non-drepped_bins_single_assemblies.fasta in1=/home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_1_trimmed.fastq in2=/home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_2_trimmed.fastq out="$element"_mapped.sam outu="$element"_unmapped_reads_combined.fastq
reformat.sh in="$element"_unmapped_reads_combined.fastq out1="$element"_unmapped_reads_R1.fq out2="$element"_unmapped_reads_R2.fq
done

#How to find total number of reads mapped:
samtools view -c -F 4 <bam>
