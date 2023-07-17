#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --time=14-00:00:00
#SBATCH --mem=150gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/read_mapping_genes_output

# output reads_per_base
coverm contig --proper-pairs-only -m reads_per_base --bam-files /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/read_mapping_genes_output/*.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt 
#output min covered 75%
coverm contig --proper-pairs-only -m mean --bam-files /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/read_mapping_genes_output/*.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0.75 --output-format dense --output-file coverm_min75.txt &> min75_stats.txt
#output counts
coverm contig --proper-pairs-only -m count --bam-files /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/gene_resolved_metaG_analyses/2500_size_filtered/read_mapping_genes_output/*.sorted.bam --min-covered-fraction 0 -t 20 --output-format dense &> coverM_counts.txt #Get counts for geTMM
echo "coverM coverage, depth, and counts finished here"
echo "finished mapping scripts"
