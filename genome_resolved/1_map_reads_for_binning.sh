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
cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/99per_sorted_bam_for_binning
# pull scaffolds >= 2500 bp
pullseq.py -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/"$element"_megahit_out/"$element"_renamed-assembly_contigs.fa -m 2500 -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/"$element"_megahit_out/"$element"_final.contigs_2500.fa
# map to these scaffolds to get sam file
bbmap.sh -Xmx48G threads=20 overwrite=t ref=/home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/"$element"_megahit_out/"$element"_final.contigs_2500.fa in1=/home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_1_trimmed.fastq in2=/home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_2_trimmed.fastq out="$element"_2500_mapped.sam
# convert sam to bam, and sort
samtools view -@ 20 -bS "$element"_2500_mapped.sam > "$element"_2500_mapped.bam
samtools sort -T "$element"_2500_mapped.sorted -o "$element"_2500_mapped.sorted.bam "$element"_2500_mapped.bam -@ 20
# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g minidfilter=0.99 in="$element"_2500_mapped.sorted.bam out="$element"_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t
#bin
runMetaBat.sh /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/"$element"_megahit_out/"$element"_final.contigs_2500.fa /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/99per_sorted_bam_for_binning/"$element"_2500_mapped99per.sorted.bam
done
