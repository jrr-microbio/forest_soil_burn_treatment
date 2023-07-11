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

### Calculate species accumuation curve, richness and shannon's H (a helpful guide: https://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html , and https://rpubs.com/an-bui/vegan-cheat-sheet)
accumcurve = specaccum(genome, method = "random", permutations=999)

plot(accumcurve, add=FALSE, ci.type="poly", col="black", lwd=2, ci=0.2, ci.lty=1, ci.col="#F8766D", ylim=c(0,100), xlab = "Samples", ylab = 'mag OTUs')

mean(specnumber(genome))
sd(specnumber(genome))
#boxplot(specnumber(genome) ~ chem$sw_pw, ylab = "# of species", col = c("#F8766D", "#00BFC4"))
#t.test(specnumber(genome) ~ chem$sw_pw)

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
anosim(Ord_dist, chem$burnfreq, permutations = 999)

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
data.scores$chemistry = chem$percentwoodyveg
#en = envfit(NMDS_Bray_genome, log_chem, permutations = 999, na.rm = TRUE)
#en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
#en_coord_cont$Compound = rownames(en_coord_cont)

###########Now do the final plot with the loadings as bars.
distance <- function(x, y, home = c(0,0)) {
  sqrt((x-home[1])^2 + (y-home[2])^2)
}

#fit_log <- envfit(NMDS_Bray_genome, log_chem, perm = 999, na.rm=TRUE)

#create loadings df for plot arrows
#v_scrs <- as.data.frame(scores(en, display = "vectors"))
#v_scrs <- cbind(v_scrs*ordiArrowMul(en), Compound = rownames(v_scrs))

#v_scrs <- v_scrs %>%
#  mutate(dist_pc1pc2 = distance(NMDS1, NMDS2), #distance in space from center
#         scaled_dist = dist_pc1pc2 * 10)

#get top 10 compounds based on distance in space from center
#top_dist_compounds <- v_scrs %>%
#  arrange(desc(scaled_dist)) %>%
#  slice_head(n = 10) %>%
#  mutate(num_lab = 1:10) %>%
#  select(Compound, num_lab)

#make plotting columns for loadings df
#v_scrs <- v_scrs %>%
#  mutate(top_c = ifelse(Compound %in% top_dist_compounds$Compound,"top","other"),
#         alpha_factor = ifelse(top_c == "top", 1, 0.5),
#         plot_pc1 = ifelse(top_c == "top", NMDS1, 0),
#        plot_pc2 = ifelse(top_c == "top", NMDS2, 0),
#         top_label = ifelse(top_c == "top", Compound, "")) %>%
#  left_join(top_dist_compounds, by = c("Compound"))


#Now plot on ggplot
plot_genomes <- data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  #geom_segment(data = rev(v_scrs))
              # ,aes(x= 0, xend = NMDS1*10, 
               #    y = 0, yend = NMDS2*10, 
                #   alpha = alpha_factor,
                 #  color = Compound),
  geom_point(aes(fill = chemistry), size = 8, shape = 21) 
  #stat_ellipse(aes(fill = chemistry, color = chemistry), alpha = 0.2, geom = "polygon", type="t") +
  #geom_mark_ellipse(aes(fill=chemistry, color = chemistry, label = chemistry))
  
  #geom_text(data = v_scrs,
   #         aes(label = num_lab,
    #            x = NMDS1*10,
     #           y = NMDS2*10),
      #      size = 3) +
  #new_scale_color()+
   #+
 # stat_ellipse(aes(fill = chemistry, color = chemistry), alpha = 0.2, geom = "polygon", type="t") +
# theme_bw() 

  plot_genomes

  
