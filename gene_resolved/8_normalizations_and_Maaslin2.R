library(tidyverse)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/gene_resolved_coverM_output/server_filtered_output/")

counts_df = read.csv("strict_mapping_table_97id_75cov_3xdepth.csv", row.names=1)

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
write.csv(relabund_from_counts, file="forest_fire_rel_abunds_97ID_75cov_3xdepth.csv")
relabund_from_counts_subset_test=head(relabund_from_counts, 200)
write.csv(relabund_from_counts_subset_test, file="forest_fire_rel_abunds_97ID_75cov_3xdepth_fakeSubset.csv")

#What im going to do now is bring in the annotation information and merge these genes into functional categories at the gene-id level. I will be using DRAM product output for this.

library(Maaslin2)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/gene_resolved_coverM_output/server_filtered_output/")
