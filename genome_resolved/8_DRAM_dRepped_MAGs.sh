#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
#source /opt/Miniconda2/miniconda2/bin/activate DRAM2


DRAM.py annotate -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/"*.fa" -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/output_DRAM_drepped_1.4.4 --min_contig_size 2500 --threads 20 --use_uniref --use_camper --use_fegenie --use_sulfur --use_vogdb &> log_trivedi_dram_drepped_1.4.4.txt

DRAM.py distill -i /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/output_DRAM_drepped_1.4.4/annotations.tsv -o /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds/distill
