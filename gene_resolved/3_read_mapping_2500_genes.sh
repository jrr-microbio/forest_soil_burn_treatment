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

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered

bowtie2-build /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/2500_genes_clustered.fna gene_DB --large-index --threads 20 #create bowtie2 index with a concatenated genomes fasta.

for element in $(<$1)
do
echo "begin bowtie2 '$element' mapping"
mkdir "$element"_mapping
cd "$element"_mapping
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 20 -x /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/gene_DB -S "$element".sam -1 /home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_1_trimmed.fastq -2 /home/projects-wrighton/Trivedi_soil_fire/reads/processed_reads/"$element"_2_trimmed.fastq #run bowtie2 with indexed genomes. This is basically preset of --fast option specified in bowtie2 manual. http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
echo "begin samtools"
samtools view -@ 20 -bS "$element".sam > "$element".bam #sam to bam
reformat.sh -Xmx100g idfilter=0.95 pairedonly=t primaryonly=t in="$element".bam out="$element"_mapped_95ID.FILTERED.bam #this will then give you the actual 95% ID mapping.
mv "$element".95ID.sorted.bam /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/read_mapping_genes_output 

echo "Finished mapping and samtools for '$element' in list" 
#clean up the files
#rm -r "$element".sam
#rm -r "$element".bam
#rm -r "$element"_mapped_95ID.FILTERED.bam
cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered #move back up to root and restart the process.
done
