library(reshape2)
library(tidyverse)
library(reticulate)
library(forcats)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/")

metabolism_matrix=read.csv("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/MAG_OTU_linking/annotations/13.1_linked_zOTU_MAGs_DRAM_wTax_presabs_v2_occupancy.csv",header=T)

#remove the pw/sw indicator status. I'll put this back later.
metabolism_matrix_short=metabolism_matrix[,-c(2:4)]
#Now convert this matrix into long format for ggplot2
metabolism_long = melt(data = metabolism_matrix_short)
#panda_metabolism_long=pandas.melt(metabolism_matrix_short, ignore_index=FALSE)

#index the factors that we want to join back in after the melt.
factors.join=metabolism_matrix %>%
  select("genome", "taxonomy", "tax.index")

metabolism.long.joined = metabolism_long %>%
  left_join(factors.join)

#r was treating the columns as numeric, so I convert them to factors since it just doesn't like numbers apparently. But this fixes that.
#metabolism.long.joined<- transform(metabolism.long.joined, genome_id= reorder(genome_id, sumpaths.index))

#metabolism.long.joined<- transform(metabolism.long.joined, genome= reorder(genome, taxonomy))

#making a color palette because that default was very hard to look at.
pokecolorpal=c("#1E77B3", "#F87F0E","#00D3FB","#FFB4EF","#2C9F2C", "#D52728","#9367BC","#8A564B", "#B2B2B2", "#E277C1", "#B4E0A8", "#BBBC23", "#32BDCE", "#6CF77A", "#FBD300", "#A0A0A0", "#8DD4C4", "#E3A795")#, "#FF7BAC", "#AAA6F0","#FF00FF")

metabolism.long.joined %>%
  filter(value > 0, 'Presence / Absence' > 1) %>%
  ggplot(aes(x=fct_reorder(genome,tax.index, .desc = TRUE), y = variable, size = value)) + 
  geom_point(aes(colour=factor(taxonomy))) + 
  scale_colour_manual(values=pokecolorpal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=1))
