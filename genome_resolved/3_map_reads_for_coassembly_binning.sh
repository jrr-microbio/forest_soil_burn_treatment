#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out

sed "s/>/>coassembly_/" /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out/final.contigs.fa > coassembly_renamed-assembly_contigs.fa

pullseq.py -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out/coassembly_renamed-assembly_contigs.fa -m 2500 -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out/coassembly_renamed-assembly_2500_contigs.fa

# map to these scaffolds to get sam file
bbmap.sh -Xmx48G threads=20 overwrite=t ref=/home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out/coassembly_renamed-assembly_2500_contigs.fa in1=/home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/all_unmapped_reads_R1.fastq in2=/home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/all_unmapped_reads_R2.fastq out=coassembly_2500_mapped.sam
# convert sam to bam, and sort
samtools view -@ 20 -bS coassembly_2500_mapped.sam > coassembly_2500_mapped.bam
samtools sort -T coassembly_2500_mapped.sorted -o coassembly_2500_mapped.sorted.bam coassembly_2500_mapped.bam -@ 20
# filter for high quality mapping (this step can be tunable and optional)
reformat.sh -Xmx100g minidfilter=0.99 in=coassembly_2500_mapped.sorted.bam out=coassembly_2500_mapped99per.sorted.bam pairedonly=t primaryonly=t
#bin
runMetaBat.sh /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out/coassembly_renamed-assembly_2500_contigs.fa /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/coassembly_unmapped_reads_to_bins/unmapped_reads_coassembly_megahit_out/coassembly_2500_mapped99per.sorted.bam
