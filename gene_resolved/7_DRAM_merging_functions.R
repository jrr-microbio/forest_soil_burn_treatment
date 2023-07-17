library(dplyr)
library(tidyr)
library(tibble)
library(janitor)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/gene_resolved_coverM_output/server_filtered_output/")

# read in annotation ids
ids= read.csv("fake_ANNOTS.csv",header=TRUE)
# create new variable "primary", selecting one id per gene. priority cazy > merops > kegg > vogdb
ids$cazy_hits_summ=as.character(ids$cazy_ids) #bring in cazy id
ids$peptidase_id_summ=as.character(ids$peptidase_id) #bring in peptidase id
ids$kegg_id_summ=as.character(ids$ko_id) #bring in ko id
ids$vogdb_hit_summ=as.character(ids$vogdb_id) #bring in vogdb id
ids$scaffold_ids=as.character(ids$scaffold) #bring in scaffold id
ids$primary=ifelse(ids$cazy_hits_summ!="",ids$cazy_hits_summ,ifelse(ids$peptidase_id_summ!="",ids$peptidase_id_summ,ifelse(ids$kegg_id_summ!="",ids$kegg_id_summ,ifelse(ids$vogdb_hit_summ!="",ids$vogdb_hit_summ,"none"))))
ids_only = ids %>% select(scaffold_ids,primary) #if else statement to give CAZY, if not MEROPS, if not KEGG, if not vOGDB, if not "none".
rm(ids)#just clear the env.
annotation_ids = ids_only #make this cleaner name
rm(ids_only) #just clear the env.

# pulled DRAM module IDs from the genome summary form in the DRAM Github.
module_info = read.table("DRAM_genome_summary_form.tsv", sep = "\t", header=T, na.strings = c("","NA"))
joined_annots_and_module=left_join(annotation_ids, module_info,by=c("primary"="gene_id")) #merge the two files.
#sum(duplicated(joined_annots_and_module$scaffold_ids)) #check no duplicates just in case

#Now we want to bring in the abundances. I will left join each gene.
abunds = read.csv("forest_fire_rel_abunds_95ID_75cov_3xdepth.csv")
#sum(duplicated(abunds$scaffold_ids)) #check no duplicates just in case
joined_abunds_annots_and_module=left_join(joined_annots_and_module, abunds, by=c("scaffold_ids"="id"))
sum(is.na(joined_abunds_annots_and_module)) #there are a bunch of NA genes. Do I just remove? Hm.
short_abunds_annot_mod <- joined_abunds_annots_and_module %>%
  select(-c(gene_description, sheet, header, subheader, potential_amg, module)) #Shorten this to just primary ID, the module, and the abunds.
short_abunds_annot_mod_no_none = short_abunds_annot_mod[!(short_abunds_annot_mod$primary == "none"),] #remove genes that could not be annotated. I'll run this both ways: With non-annotated and with only annotated.
short_abunds_noscaff_ids_for_consolidate = short_abunds_annot_mod_no_none[,-1]

pivoted_df = short_abunds_noscaff_ids_for_consolidate %>%
  pivot_longer(cols = -1,
             values_to = "abund",
             names_to="sample")

pivoted_df$NewIds <- paste(pivoted_df$primary, pivoted_df$sample, sep = "_") #for the life of me i could not figure out a simpler way to do this beasides just merging the id and the sample before the aggregate and then just splitting it up. whatever.

summed_gene_id_abunds_persamp = aggregate(abund ~ NewIds, data = pivoted_df, FUN = sum) #aggregate now if the sample/gene id are identical.
summed_gene_id_abunds_persamp <- separate(summed_gene_id_abunds_persamp, NewIds, into = c("gene_id", "sample"), sep = "_") #now just separate this again so i can make a df for Maaslin2.
Maaslin_wide_input = summed_gene_id_abunds_persamp %>% pivot_wider(names_from = "gene_id",
                                                                   values_from = c("abund")) #Pivot wider format required for Maaslin2.

Maaslin_wide_input_2 = as.data.frame(t(Maaslin_wide_input))
Maaslin_wide_input_2 = Maaslin_wide_input_2 %>%
  row_to_names(row_number = 1) #This is to make it match exactly their file "pathabundance_relab" from here: download.file("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")
row_names <- row.names(Maaslin_wide_input_2) #This needs to be numeric and for some reason isnt so im just going to fix. Extract row names first.
num_Maaslin_wide_input = as.data.frame(lapply(Maaslin_wide_input_2, as.numeric), stringsAsFactors = F) #make all numeric
row.names(num_Maaslin_wide_input) <- row_names #Now add the rownames back in.

#Test this with the Huttenhower example data.
input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
fit_data_example = Maaslin2(input_data = input_data, 
                            input_metadata = input_metadata, 
                            min_prevalence = 0,
                            normalization  = "TMM",
                            transform = "NONE",
                            output         = "demo_output_example", 
                            fixed_effects  = c("diagnosis", "dysbiosis"),
                            reference      = c("diagnosis,nonIBD"))

#Now let's do it with our own data. Import the metadata table.
forest_soil_metadata = read.table("E133_sandra_combined_sample_metadata_for_JOSUE_METAG_v2.txt", sep = "\t", header=T, row.names=1)

#And now run the same commands.
fit_data = Maaslin2(input_data = num_Maaslin_wide_input, 
                               input_metadata = forest_soil_metadata, 
                               min_prevalence = 0,
                               normalization  = "TSS",
                               transform = "LOG",
                               output         = "testing_subset_output", 
                               fixed_effects  = "burn_control")
