library(dplyr)
library(tidyr)
library(tibble)
library(janitor)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/gene_resolved_coverM_output/MaAslin2/")

# read in annotation ids
ids= read.delim("annotations_DRAM_noMinSize.tsv", header=TRUE, sep = "\t", fill = T)
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
module_info = read.delim("DRAM_genome_summary_form.tsv", sep = "\t", header=T, na.strings = c("","NA"))
joined_annots_and_module=left_join(annotation_ids, module_info,by=c("primary"="gene_id"), relationship = "many-to-many") #merge the two files.
#sum(duplicated(joined_annots_and_module$scaffold_ids)) #check no duplicates just in case

#Now we want to bring in the abundances. I will left join each gene.
abunds = read.csv("forest_fire_rel_abunds_95ID_75cov_3xdepth.csv")
#sum(duplicated(abunds$scaffold_ids)) #check no duplicates just in case
joined_abunds_annots_and_module=left_join(joined_annots_and_module, abunds, by=c("scaffold_ids"="id"))
sum(is.na(joined_abunds_annots_and_module)) #there are a bunch of NA genes. Do I just remove? Hm.
short_abunds_annot_mod <- joined_abunds_annots_and_module %>%
  select(-c(gene_description, sheet, header, subheader, potential_amg, module)) #Shorten this to just primary ID, the module, and the abunds.
short_abunds_annot_mod_no_none = short_abunds_annot_mod[!(short_abunds_annot_mod$primary == "none"),] #remove genes that could not be annotated confidently.
short_abunds_noscaff_ids_for_consolidate = short_abunds_annot_mod_no_none[,-1]
short_abunds_noscaff_ids_for_consolidate$primary <- gsub("_", "-", short_abunds_noscaff_ids_for_consolidate$primary)#i'm underscore delimiting some text so i'm just replacing underscores with dashes just in case.

pivoted_df = short_abunds_noscaff_ids_for_consolidate %>%
  pivot_longer(cols = -1,
             values_to = "abund",
             names_to="sample")

pivoted_df = pivoted_df[!(is.na(pivoted_df$primary) | pivoted_df$primary==""), ] #Removing all blank cells to be safe.

pivoted_df$NewIds <- paste(pivoted_df$primary, pivoted_df$sample, sep = "_") #for the life of me i could not figure out a simpler way to do this besides just merging the id and the sample before the aggregate and then just splitting it up. whatever.

summed_gene_id_abunds_persamp = aggregate(abund ~ NewIds, data = pivoted_df, FUN = sum) #aggregate now if the sample/gene id are identical.

summed_gene_id_abunds_persamp <- separate(summed_gene_id_abunds_persamp, NewIds, into = c("gene_id", "sample"), sep = "_") #now just separate this again so i can make a df for Maaslin2.
Maaslin_wide_input = summed_gene_id_abunds_persamp %>% pivot_wider(names_from = "gene_id",
                                                                   values_from = c("abund")) #Pivot wider format required for Maaslin2.

Maaslin_wide_input_2 = as.data.frame(t(Maaslin_wide_input))
Maaslin_wide_input_2 = Maaslin_wide_input_2 %>%
  row_to_names(row_number = 1) #Making row 1 the column to make it match exactly their file "pathabundance_relab" from here: download.file("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")
row_names <- row.names(Maaslin_wide_input_2) #This needs to be numeric and for some reason isnt so im just going to fix. Extract row names first.
num_Maaslin_wide_input = as.data.frame(lapply(Maaslin_wide_input_2, as.numeric), stringsAsFactors = F) #make all numeric
row.names(num_Maaslin_wide_input) <- row_names #Now add the rownames back in.

#Test this with the Huttenhower example data.
#input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
#input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
#fit_data_example = Maaslin2(input_data = input_data, 
#                            input_metadata = input_metadata, 
 #                           min_prevalence = 0,
  #                          normalization  = "TMM",
   #                         transform = "NONE",
    #                        output         = "demo_output_example", 
     #                       fixed_effects  = c("diagnosis", "dysbiosis"),
      #                      reference      = c("diagnosis,nonIBD"))

#Now i'm going to remove any genes from our dataset that are only present in a single sample. Using 10% occupancy (i.e., >1 sample)
num_Maaslin_wide_input_filtered <- num_Maaslin_wide_input %>%
  filter(rowSums(. != 0, na.rm = TRUE) > 1)

#Now let's do it with our own data. Import the metadata table.
forest_soil_metadata = read.table("E133_sandra_combined_sample_metadata_for_JOSUE_METAG_v2.txt", sep = "\t", header=T, row.names=1)

library(Maaslin2)

#And now run the same commands.
fit_data = Maaslin2(input_data = num_Maaslin_wide_input_filtered, 
                               input_metadata = forest_soil_metadata, 
                               min_prevalence = 0,
                               normalization  = "TSS",
                               transform = "LOG",
                               output         = "output_burncontrol_v3_filtered", 
                               fixed_effects  = "burn_control")

burncontrol_sig_results = read.delim(file="output_burncontrol_v3_filtered/significant_results.tsv")
ids= read.delim("annotations_DRAM_noMinSize.tsv", header=TRUE, sep = "\t", fill = T)
ids=ids[,-c(1:15)]
joined_annots_and_module_formerge=joined_annots_and_module[,-1] #remove first column with gene id

#this only really does the gene_id column for each gene that has it (which is a small subset).
burncontrol_sig_results=left_join(burncontrol_sig_results, joined_annots_and_module_formerge, by=c("feature"="primary")) #now left join them to get the gene descriptions. 
burncontrol_sig_results=burncontrol_sig_results[!duplicated(burncontrol_sig_results$feature),] #now remove the duplicated entries. 
burncontrol_sig_results=left_join(burncontrol_sig_results, joined_annots_and_module_formerge, by=c("feature"="primary")) #now left join them to get the gene descriptions. This does general gene_id.
burncontrol_sig_results=burncontrol_sig_results[!duplicated(burncontrol_sig_results$feature),] #now remove the duplicated entries. 
#################################### MEROPS
#################################### MEROPS
#################################### MEROPS
#Manually doing the other types of gene IDs we used to get more information.
MER0subset_burncontrol_sig_results <- subset(burncontrol_sig_results, grepl("MER", feature)) #make a subset of all merops
MER0subset_burncontrol_sig_results=left_join(MER0subset_burncontrol_sig_results, ids, by=c("feature"="peptidase_id")) #now left join them to get the gene descriptions. 
MER0subset_burncontrol_sig_results=MER0subset_burncontrol_sig_results[!duplicated(MER0subset_burncontrol_sig_results$feature),] #now remove the duplicated entries. 
#################################### VOG
#################################### VOG
#################################### VOG
VOGsubset_burncontrol_sig_results <- subset(burncontrol_sig_results, grepl("VOG", feature)) #make a subset of all VOGDB
VOGsubset_burncontrol_sig_results=left_join(VOGsubset_burncontrol_sig_results, ids, by=c("feature"="vogdb_id")) #now left join them to get the gene descriptions. 
VOGsubset_burncontrol_sig_results=VOGsubset_burncontrol_sig_results[!duplicated(VOGsubset_burncontrol_sig_results$feature),] #now remove the duplicated entries. 
#################################### KEGG
#################################### KEGG
#################################### KEGG
KEGGsubset_burncontrol_sig_results <- subset(burncontrol_sig_results, grepl("K", feature)) #make a subset of all kegg
KEGGsubset_burncontrol_sig_results=left_join(KEGGsubset_burncontrol_sig_results, ids, by=c("feature"="ko_id")) #now left join them to get the gene descriptions. 
KEGGsubset_burncontrol_sig_results=KEGGsubset_burncontrol_sig_results[!duplicated(KEGGsubset_burncontrol_sig_results$feature),] #now remove the duplicated entries. 
#################################### GH
#################################### GH
#################################### GH
GHsubset_burncontrol_sig_results <- subset(burncontrol_sig_results, grepl("GH", feature)) #make a subset of all kegg
GHsubset_burncontrol_sig_results=left_join(GHsubset_burncontrol_sig_results, ids, by=c("feature"="cazy_ids")) #now left join them to get the gene descriptions. 
GHsubset_burncontrol_sig_results=GHsubset_burncontrol_sig_results[!duplicated(GHsubset_burncontrol_sig_results$feature),] #now remove the duplicated entries. 
#Now im just going to write these files out.
write.csv(MER0subset_burncontrol_sig_results, file="MER_subset_burncontrol_sig.csv")
write.csv(VOGsubset_burncontrol_sig_results, file="VOG_subset_burncontrol_sig.csv")
write.csv(KEGGsubset_burncontrol_sig_results, file="KEGG_subset_burncontrol_sig.csv")
write.csv(GHsubset_burncontrol_sig_results, file="GH_subset_burncontrol_sig.csv")

#####################
#####################
#####################

fit_data = Maaslin2(input_data = num_Maaslin_wide_input_filtered, 
                    input_metadata = forest_soil_metadata, 
                    min_prevalence = 0,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "output_burnfreq_v3_filtered", 
                    fixed_effects  = "burnfreq")

burnfreq_sig_results = read.delim(file="output_burnfreq_v3_filtered/significant_results.tsv")
#################################### MEROPS
#################################### MEROPS
#################################### MEROPS
#Manually doing the other types of gene IDs we used to get more information.
MER0subset_burnfreq_sig_results <- subset(burnfreq_sig_results, grepl("MER", feature)) #make a subset of all merops
MER0subset_burnfreq_sig_results=left_join(MER0subset_burnfreq_sig_results, ids, by=c("feature"="peptidase_id")) #now left join them to get the gene descriptions. 
MER0subset_burnfreq_sig_results=MER0subset_burnfreq_sig_results[!duplicated(MER0subset_burnfreq_sig_results$feature),] #now remove the duplicated entries. 
#################################### VOG
#################################### VOG
#################################### VOG
VOGsubset_burnfreq_sig_results <- subset(burnfreq_sig_results, grepl("VOG", feature)) #make a subset of all VOGDB
VOGsubset_burnfreq_sig_results=left_join(VOGsubset_burnfreq_sig_results, ids, by=c("feature"="vogdb_id")) #now left join them to get the gene descriptions. 
VOGsubset_burnfreq_sig_results=VOGsubset_burnfreq_sig_results[!duplicated(VOGsubset_burnfreq_sig_results$feature),] #now remove the duplicated entries. 
#################################### KEGG
#################################### KEGG
#################################### KEGG
KEGGsubset_burnfreq_sig_results <- subset(burnfreq_sig_results, grepl("K", feature)) #make a subset of all kegg
KEGGsubset_burnfreq_sig_results=left_join(KEGGsubset_burnfreq_sig_results, ids, by=c("feature"="ko_id")) #now left join them to get the gene descriptions. 
KEGGsubset_burnfreq_sig_results=KEGGsubset_burnfreq_sig_results[!duplicated(KEGGsubset_burnfreq_sig_results$feature),] #now remove the duplicated entries. 
#################################### GH
#################################### GH
#################################### GH
GHsubset_burnfreq_sig_results <- subset(burnfreq_sig_results, grepl("GH", feature)) #make a subset of all kegg
GHsubset_burnfreq_sig_results=left_join(GHsubset_burnfreq_sig_results, ids, by=c("feature"="cazy_ids")) #now left join them to get the gene descriptions. 
GHsubset_burnfreq_sig_results=GHsubset_burnfreq_sig_results[!duplicated(GHsubset_burnfreq_sig_results$feature),] #now remove the duplicated entries. 
#Now im just going to write these files out.
write.csv(MER0subset_burnfreq_sig_results, file="MER_subset_burnfreq_sig.csv")
write.csv(VOGsubset_burnfreq_sig_results, file="VOG_subset_burnfreq_sig.csv")
write.csv(KEGGsubset_burnfreq_sig_results, file="KEGG_subset_burnfreq_sig.csv")
write.csv(GHsubset_burnfreq_sig_results, file="GH_subset_burnfreq_sig.csv")

#####################
#####################
#####################

fit_data = Maaslin2(input_data = num_Maaslin_wide_input_filtered, 
                    input_metadata = forest_soil_metadata, 
                    min_prevalence = 0,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "output_woodyveg_v3_filtered", 
                    fixed_effects  = "percentwoodyveg")

woodyveg_sig_results = read.delim(file="output_woodyveg_v3_filtered/significant_results.tsv")
#################################### MEROPS
#################################### MEROPS
#################################### MEROPS
#Manually doing the other types of gene IDs we used to get more information.
MER0subset_woodyveg_sig_results <- subset(woodyveg_sig_results, grepl("MER", feature)) #make a subset of all merops
MER0subset_woodyveg_sig_results=left_join(MER0subset_woodyveg_sig_results, ids, by=c("feature"="peptidase_id")) #now left join them to get the gene descriptions. 
MER0subset_woodyveg_sig_results=MER0subset_woodyveg_sig_results[!duplicated(MER0subset_woodyveg_sig_results$feature),] #now remove the duplicated entries. 
#################################### VOG
#################################### VOG
#################################### VOG
VOGsubset_woodyveg_sig_results <- subset(woodyveg_sig_results, grepl("VOG", feature)) #make a subset of all VOGDB
VOGsubset_woodyveg_sig_results=left_join(VOGsubset_woodyveg_sig_results, ids, by=c("feature"="vogdb_id")) #now left join them to get the gene descriptions. 
VOGsubset_woodyveg_sig_results=VOGsubset_woodyveg_sig_results[!duplicated(VOGsubset_woodyveg_sig_results$feature),] #now remove the duplicated entries. 
#################################### KEGG
#################################### KEGG
#################################### KEGG
KEGGsubset_woodyveg_sig_results <- subset(woodyveg_sig_results, grepl("K", feature)) #make a subset of all kegg
KEGGsubset_woodyveg_sig_results=left_join(KEGGsubset_woodyveg_sig_results, ids, by=c("feature"="ko_id")) #now left join them to get the gene descriptions. 
KEGGsubset_woodyveg_sig_results=KEGGsubset_woodyveg_sig_results[!duplicated(KEGGsubset_woodyveg_sig_results$feature),] #now remove the duplicated entries. 
#################################### GH
#################################### GH
#################################### GH
GHsubset_woodyveg_sig_results <- subset(woodyveg_sig_results, grepl("GH", feature)) #make a subset of all kegg
GHsubset_woodyveg_sig_results=left_join(GHsubset_woodyveg_sig_results, ids, by=c("feature"="cazy_ids")) #now left join them to get the gene descriptions. 
GHsubset_woodyveg_sig_results=GHsubset_woodyveg_sig_results[!duplicated(GHsubset_woodyveg_sig_results$feature),] #now remove the duplicated entries. 
#Now im just going to write these files out.
write.csv(MER0subset_woodyveg_sig_results, file="MER_subset_woodyveg_sig.csv")
write.csv(VOGsubset_woodyveg_sig_results, file="VOG_subset_woodyveg_sig.csv")
write.csv(KEGGsubset_woodyveg_sig_results, file="KEGG_subset_woodyveg_sig.csv")
write.csv(GHsubset_woodyveg_sig_results, file="GH_subset_woodyveg_sig.csv")
