#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

cd /home/projects-wrighton/Trivedi_soil_fire/reads/completed_assemblies/_bins/MQ_HQ_bins/dRep_v2.6.2_2MAGs/dereplicated_genomes/_renamed_scaffolds

gtdbtk classify_wf -x fa --genome_dir . --out_dir gtdb_03142023 --cpus 20
