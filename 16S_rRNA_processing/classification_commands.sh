#######################################################################
#######################################################################
#######################################################################
#######################################################################
#
#For training this classifier, i went into the GTDB website and downloaded the release v207. https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/
#I downloaded the files * bac120_ssu_reps_<release>.tar.gz * and * ar53_ssu_reps_r207.tar.gz * as well as their corresponding taxonomy.tsv files from
#(https://data.gtdb.ecogenomic.org/releases/release207/207.0/)
#From the ssu_reps files notes:
#        FASTA file of 16S rRNA gene sequences identified within the set of bacterial representative genomes. The longest 
#        identified 16S rRNA sequence is selected for each representative genomes. The assigned taxonomy reflects the 
#        GTDB classification of the genome. Sequences are identified using nhmmer with the 16S rRNA model (RF00177) from the 
#        RFAM database. Only sequences with a length >=200 bp and an E-value <= 1e-06 are reported. In a small number of cases, 
#        the 16S rRNA sequences are incongruent with this taxonomic assignment as a result of contaminating 16S rRNA sequences.
#
#From the taxonomy files note:
#    GTDB taxonomy for all bacterial genomes assigned to a GTDB species cluster.
#
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#
#I then merged both of those files into a single file for training: bac120_and_ar53_ssu_reps_r207.fna and bac120_and_ar53_taxonomy_r207.tsv
#
#Step 1: Import files into qiime2.
#
#
source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 1.1_bac120_and_ar53_ssu_reps_r207.fna \
  --output-path GTDB_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 1.2_bac120_and_ar53_taxonomy_r207.tsv \
  --output-path gtdb_ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads GTDB_otus.qza \
  --i-reference-taxonomy gtdb_ref-taxonomy.qza \
  --o-classifier gtdb_classifier.qza

  qiime feature-classifier classify-sklearn \
  --i-classifier gtdb_classifier.qza \
  --i-reads 1.3_fire_soil_E133_and_SANDRA_5_667_ALL_SAMPLES_08_zOTUs_040722_TAXONOMY_HEADERS.fa \
  --o-classification gtdb_taxonomy_leaundrareads.qza

qiime metadata tabulate \
  --m-input-file gtdb_taxonomy_leaundrareads.qza \
  --o-visualization gtdb_taxonomy.qzv
