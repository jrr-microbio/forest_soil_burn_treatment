library(ggplot2)
library(vegan)
library(grid)
library(MASS)
library(RColorBrewer)
library(purrr)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(mds)
library(ggdark)
library(tidyverse)
library(patchwork)
library(wesanderson)
library(ggforce)

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/gene_resolved/")

##read in feature table with species as columns and samples as rows
gene<-t(read.csv('7.3_and_8.1_forest_fire_rel_abunds_95ID_75cov_3xdepth.csv', sep=",", header=T,check.names=TRUE, row.names = 1, stringsAsFactors = F))
#gene = as.data.frame(lapply(gene, as.numeric), stringsAsFactors = F) #make all numeric
#Now i'm going to remove any genes from our dataset that are only present in a single sample. Using 10% occupancy (i.e., >1 sample)
gene=as.data.frame(gene)
gene=t(gene)
gene=as.data.frame(gene)
gene <- gene %>%
  filter(rowSums(. != 0, na.rm = TRUE) > 1) #remove anything that is present in only 1 sample.

gene_unburn = as.data.frame(t(gene[6:9,]))
gene_burn = as.data.frame(t(gene[-c(6:9),]))

gene_burnfreq0 = as.data.frame(t(gene[6:9,]))
gene_burnfreq05 = as.data.frame(t(gene[c(4,5,10,11),]))
gene_burnfreq08 = as.data.frame(t(gene[1:3,]))

##read in chemistry (nona)
chem = read.table('7.4_and_8.2_E133_sandra_combined_sample_metadata_for_JOSUE_METAG_v2.txt', sep = '\t', header = TRUE, check.names = T)
chem = chem[-1,]
#chem = chem[-9,]
rownames(chem)=chem[,1]
chem = chem[,-1]
chem = chem[,-1]
chem$burnunit = as.factor(chem$burnunit) 
chem$burnfreq = as.factor(chem$burnfreq) 
chem$yearssinceburn = as.factor(chem$yearssinceburn) 
chem$treatment = as.factor(chem$treatment) 

###############################################################################
# Richness and shannon's stats
###############################################################################
mean(specnumber(gene))
sd(specnumber(gene))
mean(specnumber(t(gene_unburn)))
sd(specnumber(t(gene_unburn)))
mean(specnumber(t(gene_burn)))
sd(specnumber(t(gene_burn)))
boxplot(specnumber(gene) ~ chem$treatment, ylab = "# of species", col = c("#F8766D", "#00BFC4"))
t.test(specnumber(gene) ~ chem$treatment)

boxplot(specnumber(gene) ~ chem$burnfreq, ylab = "# of species", col = c("#c5c5c5", "#ff721d", "#b0091a"))
boxplot(specnumber(gene) ~ chem$yearssinceburn, ylab = "# of species", col = c("#c5c5c5", "#8ad2db"))
t.test(specnumber(gene) ~ chem$yearssinceburn)


mean(diversity(t(gene_unburn), index = "shannon"))
mean(diversity(t(gene_burn), index = "shannon"))
sd(diversity(t(gene_unburn), index = "shannon"))
sd(diversity(t(gene_burn), index = "shannon"))
boxplot(diversity(gene) ~ chem$treatment, ylab = "Shannon's H'", col = c("#F8766D", "#00BFC4")) 
t.test(diversity(gene) ~ chem$treatment)

boxplot(diversity(gene) ~ chem$burnfreq, ylab = "Shannon's H'", col = c("#c5c5c5", "#ff721d", "#b0091a"))
boxplot(diversity(gene) ~ chem$yearssinceburn, ylab = "Shannon's H'", col = c("#c5c5c5", "#8ad2db"))
t.test(diversity(gene) ~ chem$yearssinceburn)

###############################################################################
# NMDSr - burn vs unburn
###############################################################################

#This establishes the coordinates for the samples and all OTUS onto an non-metric dimensional scaling
Ord_dist <-metaMDSdist(gene, distance = "bray", noshare = 0.1, trace = 1, autotransform=T)

#Run mrpp and ANOSIM treatments between samples
mrpp(gene, chem$burnunit, permutations=999, distance="bray")
mrpp(gene, chem$treatment, permutations=999, distance="bray")
mrpp(gene, chem$burnfreq, permutations=999, distance="bray")

#########################
##Now to work on plotting the NMDS on ggplot in 2D (k = 3 gives 3D). Next 5 lines just making sure the plot still looks good with added features on metaMDS command.
NMDS_Bray_gene <-metaMDS(gene, distance = "bray", k =2,
                           noshare = 0.1, trace = 1, trymax = 500, autotransform=T)
ord.gene = as.data.frame(scores((NMDS_Bray_gene), display="sites"))
#ord.gene$sampleid=row.names(ord.gene)

stressplot(NMDS_Bray_gene)

#Now to start actually plotting this on ggplot.
#Pull out the scores with the new R framework
data.scores = as.data.frame(scores(NMDS_Bray_gene)$sites)
data.scores$chemistry = chem$burnunit

###########Now do the final plot with the loadings as bars.
distance <- function(x, y, home = c(0,0)) {
  sqrt((x-home[1])^2 + (y-home[2])^2)
}

#Now plot on ggplot
plot_genes <- data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  #geom_segment(data = rev(v_scrs))
  # ,aes(x= 0, xend = NMDS1*10, 
  #    y = 0, yend = NMDS2*10, 
  #   alpha = alpha_factor,
  #  color = Compound),
  geom_point(aes(fill = chemistry), size = 8, shape = 21) +
  geom_mark_ellipse(aes(fill=chemistry, color = chemistry, label = chemistry)) 

plot_genes
