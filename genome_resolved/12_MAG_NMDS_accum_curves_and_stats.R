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

setwd("/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/")

##read in feature table with species as columns and samples as rows
genome<-t(read.table('/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/read_mapping/59_MAGs_3x_depth_75AF_TMM_normalized_nozero.tsv', sep="\t", header=T,check.names=TRUE, row.names = 1))

genome_t = t(genome)

genome_unburn = as.data.frame(t(genome[6:9,]))
genome_burn = as.data.frame(t(genome[-c(6:9),]))

genome_burnfreq0 = as.data.frame(t(genome[6:9,]))
genome_burnfreq05 = as.data.frame(t(genome[c(4,5,10,11),]))
genome_burnfreq08 = as.data.frame(t(genome[1:3,]))

##read in chemistry (nona)
chem = read.table('/Volumes/Macintosh HD/Users/josue.rodriguez/Library/CloudStorage/GoogleDrive-jarora2213@gmail.com/My Drive/University/wrighton_lab_phd/Trivedi_collaboration/E133_sandra_combined_sample_metadata_for_JOSUE_METAG_v2.txt', sep = '\t', header = TRUE, check.names = T)
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
# Species accumulation curves
###############################################################################

### Calculate species accumuation curve, richness and shannon's H (a helpful guide: https://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html , and https://rpubs.com/an-bui/vegan-cheat-sheet)
accumcurve = specaccum(genome, method = "random", permutations=999)
accumcurve_unburn = specaccum(t(genome_unburn), method = "random", permutations=999)
accumcurve_burn = specaccum(t(genome_burn), method = "random", permutations=999)

plot(accumcurve, add=F, ci.type="poly", col="black", lwd=2, ci=0.2, ci.lty=1, ci.col="#238b45", ylim=c(0,60), xlab = "Samples", ylab = 'mag OTUs')

plot(accumcurve_burn, add=T, ci.type="poly", col="black", lwd=2, ci=0.2, ci.lty=1, ci.col="#F8766D", ylim=c(0,45), xlab = "Samples", ylab = 'mag OTUs')

plot(accumcurve_unburn, add=T, ci.type="poly", col="black", lwd=2, ci=0.2, ci.lty=1, ci.col="#00BFC4", ylim=c(0,45), xlab = "Samples", ylab = 'mag OTUs')

###############################################################################
# Richness and shannon's stats
###############################################################################
mean(specnumber(genome))
sd(specnumber(genome))
mean(specnumber(t(genome_unburn)))
sd(specnumber(t(genome_unburn)))
mean(specnumber(t(genome_burn)))
sd(specnumber(t(genome_burn)))
boxplot(specnumber(genome) ~ chem$treatment, ylab = "# of species", col = c("#F8766D", "#00BFC4"))
t.test(specnumber(genome) ~ chem$treatment)

boxplot(specnumber(genome) ~ chem$burnfreq, ylab = "# of species", col = c("#c5c5c5", "#ff721d", "#b0091a"))
boxplot(specnumber(genome) ~ chem$yearssinceburn, ylab = "# of species", col = c("#c5c5c5", "#8ad2db"))
t.test(specnumber(genome) ~ chem$yearssinceburn)


mean(diversity(t(genome_unburn), index = "shannon"))
mean(diversity(t(genome_burn), index = "shannon"))
sd(diversity(t(genome_unburn), index = "shannon"))
sd(diversity(t(genome_burn), index = "shannon"))
boxplot(diversity(genome) ~ chem$treatment, ylab = "Shannon's H'", col = c("#F8766D", "#00BFC4")) 
t.test(diversity(genome) ~ chem$treatment)

boxplot(diversity(genome) ~ chem$burnfreq, ylab = "Shannon's H'", col = c("#c5c5c5", "#ff721d", "#b0091a"))
boxplot(diversity(genome) ~ chem$yearssinceburn, ylab = "Shannon's H'", col = c("#c5c5c5", "#8ad2db"))
t.test(diversity(genome) ~ chem$yearssinceburn)

###############################################################################
# NMDSr - burn vs unburn
###############################################################################

####Simple Ord with env factors
#chem.pca_log<-prcomp(na.omit(chem), center=TRUE, scale.=TRUE)
ordMAG<-metaMDS(genome, distance = "bray")
#fit_logMAG <- envfit(ordMAG, log_chem, perm = 999, na.rm=TRUE)
#scores(fit_logMAG, "vectors")

plot(ordMAG)

#This establishes the coordinates for the samples and all OTUS onto an non-metric dimensional scaling
Ord_dist <-metaMDSdist(genome, distance = "bray", noshare = 0.1, trace = 1, autotransform=T)

#Run mrpp and ANOSIM on sw vs pw between samples
mrpp(genome, chem$burnunit, permutations=999, distance="bray")
mrpp(genome, chem$treatment, permutations=999, distance="bray")
mrpp(genome, chem$burnfreq, permutations=999, distance="bray")

#########################
##Now to work on plotting the NMDS on ggplot in 2D (k = 3 gives 3D). Next 5 lines just making sure the plot still looks good with added features on metaMDS command.
NMDS_Bray_genome <-metaMDS(genome, distance = "bray", k =2,
                          noshare = 0.1, trace = 1, trymax = 500, autotransform=T)
ord.genome = as.data.frame(scores((NMDS_Bray_genome), display="sites"))
#ord.genome$sampleid=row.names(ord.genome)

stressplot(NMDS_Bray_genome)

#Now to start actually plotting this on ggplot.
#Pull out the scores with the new R framework
NMDS_Bray_genome <-metaMDS(genome, distance = "bray", k =2,
                              noshare = 0.1, trace = 1, trymax = 500, autotransform=T)

data.scores = as.data.frame(scores(NMDS_Bray_genome)$sites)
data.scores$chemistry = chem$treatment

###########Now do the final plot with the loadings as bars.
distance <- function(x, y, home = c(0,0)) {
  sqrt((x-home[1])^2 + (y-home[2])^2)
}

#Now plot on ggplot
plot_genomes <- data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  #geom_segment(data = rev(v_scrs))
              # ,aes(x= 0, xend = NMDS1*10, 
               #    y = 0, yend = NMDS2*10, 
                #   alpha = alpha_factor,
                 #  color = Compound),
  geom_point(aes(fill = chemistry), size = 8, shape = 21) +
  geom_mark_ellipse(aes(fill=chemistry, color = chemistry, label = chemistry)) 

  plot_genomes
  
  ###############################################################################
  # Next up - indicator - burn vs unburn
  ###############################################################################
  
  library(indicspecies)
  
  set.seed(1234) #set a random seed
  indic.burntreat <- multipatt(genome, chem$burnunit, duleg = F, func = "r.g",
                                permutations = 9999) #test here whether the genome is an indocator for whatever treatment
  
  indic.burntreat <- indic.burntreat$sign #export results from object of multipatt
  indic.burntreat$MAG_id <- rownames(indic.burntreat) #add in a column with MAG id
  indic.burntreat$p.fdr <- p.adjust(indic.burntreat$p.value, method = "fdr") #run a false discovery rate correction 
  indic.burntreat <- indic.burntreat %>%
    filter(p.fdr < 0.05)
  table(indic.burntreat$index) #filter out everything that isn't significant
