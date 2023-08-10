#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
#source /opt/Miniconda2/miniconda2/bin/activate DRAM2

#MAGs present in both treatments
DRAM.py annotate -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/2500_genes_clustered.fna -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/DRAM_clustered_genes --min_contig_size 2500 --threads 20 --use_camper --use_fegenie --use_sulfur --use_vogdb &> log_trivedi_dram_clustered_genes.txt

DRAM.py distill -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/DRAM_clustered_genes/annotations.tsv -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/DRAM_clustered_genes/distill
