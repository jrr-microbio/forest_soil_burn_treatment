#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

# pre-processing:
# concatenate MAGs into single file for mapping reference.
#Select medium-high quality genomes
# list.txt is a list of trimmed read names.

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/

bowtie2-build /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/all_MQHQ_drepped_99id_MAGs.fasta mag_DB --large-index --threads 20 #create bowtie2 index with a concatenated genomes fasta.

for element in $(<$1)
do
echo "begin bowtie2 '$element' mapping"
mkdir "$element"_mapping
cd "$element"_mapping
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 20 -x /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/mag_DB -S "$element".sam -1 /home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_1_trimmed.fastq -2 /home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_2_trimmed.fastq #run bowtie2 with indexed genomes. This is basically preset of --fast option specified in bowtie2 manual. http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
echo "begin samtools"
samtools view -@ 20 -bS "$element".sam > "$element".bam #sam to bam
reformat.sh -Xmx100g idfilter=0.97 pairedonly=t primaryonly=t in="$element".bam out="$element"_mapped_97ID.FILTERED.bam #this will then give you the actual 97% ID mapping. This was wrong before on bbmap as the 97% ID cutoffs were not actually cutoffs.
samtools sort -T "$element".97ID.sorted -o "$element".97ID.sorted.bam "$element"_mapped_97ID.FILTERED.bam -@ 20 #sort the bam file
mv "$element".97ID.sorted.bam /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/sorted_BAM_97ID #now we move out the sorted BAM into its own directory where all will be for coverM. 
echo "Finished mapping and samtools for '$element' in list" 
cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/ #move back up to root and restart the process.
done

#now that we've moved over all the sorted SAMs into their own directory, we run coverM.

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/sorted_BAM_97ID

echo "begin coverM"
# output reads_per_base
coverm genome --proper-pairs-only -m reads_per_base --genome-fasta-extension fa --genome-fasta-directory /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/ --bam-files *.sorted.bam --threads 20 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0 --output-format dense --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt 
#output min covered 75%
coverm genome --proper-pairs-only -m mean --genome-fasta-extension fa --genome-fasta-directory /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/ --bam-files *.sorted.bam --threads 20 --min-read-percent-identity-pair 0.97 --min-covered-fraction 0.75 --output-format dense --output-file coverm_min75.txt &> min75_stats.txt
#output counts
coverm genome --proper-pairs-only -m count --genome-fasta-extension fa --genome-fasta-directory /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/ --bam-files *.sorted.bam --min-covered-fraction 0 -t 20 --output-format dense &> coverM_counts.txt #Get counts for geTMM
echo "coverM coverage, depth, and counts finished here"
echo "finished mapping scripts"
