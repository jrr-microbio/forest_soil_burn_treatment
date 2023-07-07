#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --job-name=join_asv
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/qiime2_gtdb_classifier/join_asv_bins_script

#Following instructions here: https://github.com/WrightonLabCSU/join_asvbins

join_asvbins \
        -b /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/_80_dRepped_MAGs \
        -a /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/qiime2_gtdb_classifier/classifying_zOTU_sequences/fire_soil_E133_and_SANDRA_5_667_ALL_SAMPLES_08_zOTUs_040722_TAXONOMY_HEADERS.fa \
        -o /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/qiime2_gtdb_classifier/join_asv_bins_script \
        --generic_16S /home/projects-wrighton/Trivedi_soil_fire/reads/all_bins_renamed_scaffolds/qiime2_gtdb_classifier/classifying_zOTU_sequences/bac120_and_ar53_ssu_reps_r207.fna \
        -t 10
