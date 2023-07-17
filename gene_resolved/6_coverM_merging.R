#Note: in order to read the coverM files in you need to slightly modify the top rows of the file. R doesn't like them. If file is too big to open in excel, in bash you can do: 
#tail -n +15 input_file.txt > output_file.txt 
#You have to change "15" to however many lines you have. This is N+15, so it will remove the top 14 rows. If you have more or less samples, you need to modify accordingly.
library(tidyverse)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/gene_resolved_coverM_output/server_filtered_output/")
###For running R on server: /home/opt/R_4.1.2/R-4.1.2/bin/R
#Read in counts and modify the headers.
counts = read.table("95_ID_coverM_counts_notoprow.txt", sep = "\t", row.names=1)
counts[1, ] = str_replace(counts[1, ], ".95ID.sorted Read Count", "_readcounts")
#function to make the first row the header names
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
#apply the function
counts_headerfix = header.true(counts)
#reduce size for testing
#counts_reduced = counts_headerfix[1:100, ]
#Going to remove everything that is just zero across the board to reduce compute time.
#counts_headerfix$sum_counts = rowSums(counts_headerfix) #sum rows
#counts_headerfix = counts_df[!(counts_headerfix$sum_counts == 0),] #remove those equal to 0 in new counts_df$sum_counts
#counts_headerfix = subset(counts_headerfix, select = -c(counts_headerfix)) #remove the sum column now that we're done with it.

###
coverage = read.table("95_ID_coverm_min75.txt", sep = "\t", row.names=1)
coverage[1, ] = str_replace(coverage[1, ], ".95ID.sorted Mean", "_coverage")
coverage_headerfix = header.true(coverage)
#coverage_reduced = coverage_headerfix[1:100, ]

depth = read.table("95_ID_coverm_reads_per_base.txt", sep = "\t", row.names=1)
depth[1, ] = str_replace(depth[1, ], ".95ID.sorted Reads per base", "_depth")
depth_headerfix = header.true(depth)
#depth_reduced = depth_headerfix[1:100, ]

#I removed the first column because as we can see the sample was not good and it has no good hits at the specified coverage.
counts_new = counts_headerfix[, -1]
coverage_new = coverage_headerfix[, -1]
depth_new = depth_headerfix[, -1]

#counts_reduced_new = counts_reduced[, -1]
#coverage_reduced_new = coverage_reduced[, -1]
#depth_reduced_new = depth_reduced[, -1]
#Now i have the 11 samples that were good and recruited reads. I want to "merge" these. I want to keep hits that were >97% ID, >=3x depth, and >=75% coverage. Taking off kai and kayla's script:


## ---------------------------
## Script merge_coverM_outputs
## 
## Purpose merge MetaG coverM outputs
##
## By Mikayla Borton 
## Date adapted: May 1, 2023
##
## ---------------------------
## Script adapted from place_coassemblies_to_treatment script
## Author: Ikaia Leleiwi
##
## Date Created: May 16, 2022
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
## ---------------------------

##Libraries##
#library(edgeR)
library(tibble)
#mapping files
#default coverM output with min_covered_fraction set
covered_fraction <- coverage_new
#coverM output with -m reads_per_base
reads_per_base <- depth_new
#coverM output with -m trimmed_mean
trimmed_mean <- counts_new

covered_fraction_v2 = rownames_to_column(covered_fraction, "contigs")
#covered_fraction_v2[1, ] = as.character(colnames(covered_fraction_v2))
#names(covered_fraction_v2) = NULL
reads_per_base_v2 = rownames_to_column(reads_per_base, "contigs")
trimmed_mean_v2 = rownames_to_column(trimmed_mean, "contigs")

#modify longer function
#note your sample names should not have underscores, this is looking for the first underscore and then will split into sample id and map_type
mod_long <- function(df){
  df <- df %>%
    pivot_longer(cols = -1,
                 values_to = "count",
                 names_to = "sample") %>%
    separate(col = "sample",
             into = c("sample_id", "map_type"),
             extra = "merge")
}

##Combine data##
df_list <- list(covered_fraction_v2,
                trimmed_mean_v2,
                reads_per_base_v2)


combined_df <- lapply(df_list, mod_long) %>%
  reduce(rbind)

combined_df_wide <- combined_df %>%
  pivot_wider(names_from = "map_type",
              values_from = c("count"))

##Filter data##
bin_counts_matrix <- combined_df_wide %>%
  mutate(coverage = ifelse(coverage > 0, readcounts, 0), #calculate genomes with 75% MAG coverage. If cov >0, write counts, otherwise give 0
         coverage = ifelse((as.numeric(depth))*151 >= 3, readcounts, 0)) %>% #calculate 3x coverage per base. If depth*151 >=3, give counts, otherwise give 0.
  select(contigs, sample_id, coverage) %>%
  pivot_wider(values_from = coverage,
              names_from = sample_id) %>%
  arrange(contigs) %>%
  column_to_rownames(var = "contigs") %>%
  as.matrix()

bin_counts_out <- bin_counts_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "bin.id")

write_csv(bin_counts_out, "strict_mapping_table_95id_75cov_3xdepth.csv")

###################
#this is the last step that i had to do on server. Moving this over to local machine now.
###################

counts_df = read.csv("strict_mapping_table_95id_75cov_3xdepth.csv", row.names=1)

#Going to remove everything that is just zero across the board to reduce compute time.
counts_df$sum_counts = rowSums(counts_df) #sum rows
counts_df_filtered = counts_df[!(counts_df$sum_counts == 0),] #remove those equal to 0 in new counts_df$sum_counts
counts_df_filtered = subset(counts_df_filtered, select = -c(sum_counts)) #remove the sum column now that we're done with it.

#Function for calculating relative abundance
relabund <- function(df, columns = c(NA)) 
  #takes feature table and calculates relative abundance of columns, omitting NA's
  #considers only columns listed in columns argument (character vector of column names). 
  #columns (default) = all columns
{
  if (NA %in% columns){
    df <- sweep(df, 2, colSums(df, na.rm = TRUE), '/')
  }
  else {
    df_relabund <- df %>% select(all_of(columns)) %>%
      sweep(., 2, colSums(., na.rm = TRUE), '/')
    df[,columns] <- df_relabund
  }
  
  return(df)
}

relabund_from_counts = relabund(counts_df_filtered) #apply rel_abund function
write.csv(relabund_from_counts, file="forest_fire_rel_abunds_95ID_75cov_3xdepth.csv")
