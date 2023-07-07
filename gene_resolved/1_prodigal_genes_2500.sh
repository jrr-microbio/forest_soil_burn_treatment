#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=prodigal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered

#Prodigal V2.6.2: February, 2015
prodigal -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/all_assembled_scaffolds_2500.fasta -o 2500_genes.gff -a 2500_genes.faa -d 2500_genes.fna -s 2500_genes.stats -f gff -p meta
