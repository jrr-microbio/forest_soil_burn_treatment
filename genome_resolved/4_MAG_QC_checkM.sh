#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/_all_bins
checkm lineage_wf -t 20 -x fa . checkm #run checkM
cd checkm
checkm qa -o 2 -f results.txt --tab_table -t 20 lineage.ms . #run checkM table create
awk -F "\t" '{if ($6 >49 && $7 <11) print $1}' results.txt #pull out MQ / HQ bins
