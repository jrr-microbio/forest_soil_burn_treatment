#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=100gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads
mkdir processed_reads
mkdir fastqc_files

for element in $(<$1)
do 
mkdir /home/projects-wrighton/Trivedi_soil_fire/reads/fastqc_files/"$element"/
mkdir /home/projects-wrighton/Trivedi_soil_fire/reads/fastqc_files/"$element"/fastqc_pretrim/
fastqc "$element"_1.fq "$element"_2.fq -o /home/projects-wrighton/Trivedi_soil_fire/reads/fastqc_files/"$element"/fastqc_pretrim/
sickle pe -f "$element"_1.fq -r "$element"_2.fq -t sanger -o "$element"_1_trimmed.fastq -p "$element"_2_trimmed.fastq -s "$element"_R1R2_singles_trimmedout.fastq
rm "$element"_R1R2_singles_trimmedout.fastq
mkdir /home/projects-wrighton/Trivedi_soil_fire/reads/fastqc_files/"$element"/fastqc_posttrim/
fastqc "$element"_1_trimmed.fastq "$element"_2_trimmed.fastq -o /home/projects-wrighton/Trivedi_soil_fire/reads/fastqc_files/"$element"/fastqc_posttrim/
done

mv *.fastq /home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/

##we could also use this other and compare which is best- but not necessary for this dataset:
#bbduk.sh threads=10 overwrite=t in1=../raw_reads/WCRC_304_R1.fastq in2=../raw_reads/WCRC_304_R2.fastq ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 minlength=75 ref=/opt/bbtools/bbmap/resources/adapters.fa out1=WCRC_304_R1_bbduktrimmed.fastq out2=WCRC_304_R2_bbduktrimmed.fastq
