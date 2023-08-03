#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --job-name=cd-hit
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/DRAM_unclustered_genes/

#====== CD-HIT version 4.6 (built on Nov  6 2014) ======
cd-hit-est -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/DRAM_unclustered_genes/working_dir/all_assembled_scaffolds_2500/genes.annotated.fna -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/DRAM_unclustered_genes/2500_genes_clustered.fna -c 0.95 -aS 0.80 -M 128000 -T 20 -d 0 -B 1 -bak 1