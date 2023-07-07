#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
#source /opt/Miniconda2/miniconda2/bin/activate DRAM2

#MAGs present in both treatments
DRAM.py annotate -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/both/"*.fa" -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/both/output_DRAM_drepped_1.4.4_both --min_contig_size 2500 --threads 20 --use_uniref --use_camper --use_fegenie --use_sulfur --use_vogdb --gtdb_taxonomy /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/gtdb_03142023/gtdb_summary_bac-ar.tsv &> log_trivedi_dram_drepped_1.4.4_both.txt

DRAM.py distill -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/both/output_DRAM_drepped_1.4.4_both/annotations.tsv -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/both/output_DRAM_drepped_1.4.4_both/distill

#MAGs present in burned treatments
DRAM.py annotate -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/burned/"*.fa" -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/burned/output_DRAM_drepped_1.4.4_burned --min_contig_size 2500 --threads 20 --use_uniref --use_camper --use_fegenie --use_sulfur --use_vogdb --gtdb_taxonomy /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/gtdb_03142023/gtdb_summary_bac-ar.tsv &> log_trivedi_dram_drepped_1.4.4_burned.txt

DRAM.py distill -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/burned/output_DRAM_drepped_1.4.4_burned/annotations.tsv -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/burned/output_DRAM_drepped_1.4.4_burned/distill

#MAGs present in control treatments
DRAM.py annotate -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/unburned/"*.fa" -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/unburned/output_DRAM_drepped_1.4.4_unburned --min_contig_size 2500 --threads 20 --use_uniref --use_camper --use_fegenie --use_sulfur --use_vogdb --gtdb_taxonomy /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/gtdb_03142023/gtdb_summary_bac-ar.tsv &> log_trivedi_dram_drepped_1.4.4_unburned.txt

DRAM.py distill -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/unburned/output_DRAM_drepped_1.4.4_unburned/annotations.tsv -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/unburned/output_DRAM_drepped_1.4.4_unburned/distill

#MAGs present in neither treatments
DRAM.py annotate -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/neither/"*.fa" -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/neither/output_DRAM_drepped_1.4.4_neither --min_contig_size 2500 --threads 20 --use_uniref --use_camper --use_fegenie --use_sulfur --use_vogdb --gtdb_taxonomy /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/gtdb_03142023/gtdb_summary_bac-ar.tsv &> log_trivedi_dram_drepped_1.4.4_neither.txt

DRAM.py distill -i /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/neither/output_DRAM_drepped_1.4.4_neither/annotations.tsv -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/output_DRAM_drepped_burn_vs_unburned/neither/output_DRAM_drepped_1.4.4_neither/distill
