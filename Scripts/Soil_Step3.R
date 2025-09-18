################################################################################################################
### Project: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird
### Script: Soil_Step1.R
### This script is step 3 of 6
### Script objective(s): Examine whether soil chemical profiles vary across the landscape
### Script objective(s) continued: Create multiple files for export to PRIMER. Make figures from the results of the PRIMER analysis
################################################################################################################

######### Setup #############

# load packages
library(here)
library(tidyverse)
library(broom)
library(vegan)
library(perm)
library(RColorBrewer)
library(colorspace)
library(ecodist)
library(writexl)
library(geodist) 
library(patchwork)

########## Import burrow coordinates and make pairwise distance matrix #########

# import burrow coordinates
burrowcoords <- read.csv(here("Data", "burrow_coordinates.csv"), header=TRUE)
head(burrowcoords)

# make data frame with only latitude and longitude
coords <- burrowcoords %>%
  arrange(Burrow) %>% # arrange by Burrow to ensure organization of final matrix is the same as the soil chemical matrices
  dplyr::select(longitude, latitude)

# get matrix of distances between burrow coordinates
burrows <- geodist(coords) # default measure is "cheap" which is okay for small areas (size of city)
rownames(burrows) <- burrowcoords$Burrow
colnames(burrows) <- burrowcoords$Burrow
class(burrows) # matrix/array

burrows_df <- as.data.frame(burrows) # convert to data frame
# export copy for use in PRIMER
writexl::write_xlsx(burrows_df, here("PRIMER/Imported Files", "burrowdistancematrix.xlsx"))

# convert to a dist object with only the lower triangle of the matrix
burrdist <- as.dist(burrows)

hist(burrdist)
min(burrdist)

###############################################################################################################
###################  Perform a mantel test using background and burrow samples ###############################

# objective: do chemical profiles vary with geographic distance across the colony for a) background and b) burrow samples?

# prepare data for background and burrow samples
# we need the average StdPA for each compound (RT) across the 3 replicate samples (so one measure per location)

# import standardize peak areas for soil data
StdPA_Soil <- readRDS(here("Outputs", "StdPA_Soil.rds"))

# find average of each RT (compound) across the triplicate samples
AveStdPA_Soil <- StdPA_Soil %>% 
  group_by(Location, RT) %>% 
  summarize(Ave_StdPA=mean(Std_PA))

# convert to wide format
AveStdPA_Soil_wide <- pivot_wider(AveStdPA_Soil, names_from = RT, values_from= Ave_StdPA)

StdPA_Soil_id <- StdPA_Soil %>%
  dplyr::select(Location, Burrow, Type, Site) %>%
  distinct()

nrow(StdPA_Soil_id)
nrow(AveStdPA_Soil_wide)

# combine id info with wide format
AveStdPA_Soil_id <- inner_join(StdPA_Soil_id, AveStdPA_Soil_wide)

# simplify to only bkgd samples
AveStdPA_bkgd <- AveStdPA_Soil_id %>%
  filter(Type == "Bkgd") %>%
  arrange(Burrow)

# simplify to only burrow samples
AveStdPA_burr <- AveStdPA_Soil_id %>%
  filter(Type == "Burr") %>%
  arrange(Burrow)

# make a version with all samples for use in PRIMER
AveStdPA_all <- AveStdPA_Soil_id %>%
  arrange(Burrow)

# export copies for use in PRIMER
AveStdPA_burr %>% select(-Location, -Type, -Site) %>% writexl::write_xlsx(here("PRIMER/Imported Files", "AveStdPA_burr.xlsx"))
AveStdPA_bkgd %>% select(-Location, -Type, -Site) %>% writexl::write_xlsx(here("PRIMER/Imported Files", "AveStdPA_bkgd.xlsx"))
AveStdPA_all %>% select(-Location, -Type, -Site) %>% writexl::write_xlsx(here("PRIMER/Imported Files", "AveStdPA_All.xlsx"))


### mantel test for background samples

# make BC dissimilarity matrix
bkgd_ave <- AveStdPA_bkgd %>% dplyr::select(-Location, -Burrow, -Type, -Site)
bkgd_ave_bray <- vegdist(bkgd_ave[,], method="bray")

# run mantel test with 9999 permutations of the data
set.seed(682)
bkgd_mantel <-vegan::mantel( xdis= bkgd_ave_bray,ydis= burrdist, method="spearman", permutations=9999)
bkgd_mantel # view results
saveRDS(bkgd_mantel, here("Models", "background_mantel.rds")) # save model

### mantel test for burrow samples

# make BC dissimilarity matrix
burr_ave <- AveStdPA_burr %>% dplyr::select(-Location, -Burrow, -Type, -Site)
burr_ave_bray <- vegdist(burr_ave[,], method="bray")

# run mantel test with 9999 permutations of the data
set.seed(409)
burr_mantel <- vegan::mantel(xdis= burr_ave_bray, ydis= burrdist, method="spearman", permutations=9999)
burr_mantel # view results
saveRDS(burr_mantel, here("Models", "burrow_mantel.rds")) # save model

# Next, we can test whether the correlation detected using the mantel test can be strengthened using just a subset of compounds in the chemical profile
# To do this, we will use BVStep
# this method allow us to find the best subset of environmental variables that explains correlation between two matrices
# this step will be performed in PRIMER using the excel files exported above
# for the results of this analysis see: 
# PRIMER/Exported Results/BVSTEP Bkgd Compounds.rtf and PRIMER/Exported Results/BVSTEP Burrow Compounds.rtf


# perform a final Mantel test to examine the degree of covariance between the temporal and spatial aspects of the soil sampling

# get the sampling dates
dates <- burrowcoords %>%
  arrange(Burrow) %>% # arrange by Burrow to ensure organization of final matrix is the same as the soil chemical matrices
  column_to_rownames(var="Burrow") %>%
  dplyr::select(Sampling_Day)

datesdist <- dist(dates) # make a pairwise distance matrix using the sampling dates

# run a Mantel test using the spatial burrow distance matrix (burrdist) and the temporal distance matrix
set.seed(567)
date_space <- vegan::mantel(xdis= datesdist, ydis= burrdist, permutations=9999)

###############################################################################################################
##### PERMANOVA and CAP to explore differences in chemical profiles based on colony site and sample type ######

# We can also explore whether soil sample chemical profiles vary by site? or type/class? 
# site refers to geographic location (sites 1, 2, and 3)
# class refers to background, occupied burrow, and unoccupied burrow
# we will examine this using PERMANOVA and CAP (Canonical analysis of principal coordinates)
# these are unconstrained and constrained multivariate analyses, respectively 
# we have an unbalanced design due to not having equal numbers of samples per Site or Class
# the methods available in R are not optimized to handle unbalanced designs
# we will export data to PRIMER and use CAP and PERMANOVA methods available there

# prep data for export
StdPA_Soil_permanova <- StdPA_Soil %>%
  mutate(Class = ifelse(Type == "Bkgd", "Bkgd",
                        ifelse(Type== "Burr" & Occupied == "Y", "Occ_Burrow", "Unocc_Burrow")))%>%
  dplyr::select(Location, Burrow, Site, Class) %>%
  distinct()

permanova <- inner_join(StdPA_Soil_permanova, AveStdPA_Soil_wide, by="Location")
writexl::write_xlsx(permanova, here("PRIMER/Imported Files", "AveStdPA_All.xlsx"))

# see "PRIMER/Exported Results/Permanova Site and Class.rft" for PERMANOVA results exported from PRIMER
# see "PRIMER/Exported Results/CAP for Class.rft" for CAP results testing importance of Class for differentiating samples
# see "PRIMER/Exported Results/CAP for Site.rft" for CAP results testing importance of Site for differentiating samples

# make figures to go along with the CAP analysis
# to do this, we will import scores (aka coordinates) assigned to the samples by the CAP analyses
# these were exported from PRIMER
CAPscores_class <- read.csv(here("PRIMER/Exported Results", "CAP Scores Class.csv"), header=T) # for landscape features
CAPscores_site <- read.csv(here("PRIMER/Exported Results", "CAP Scores Site.csv"), header=T) # for colony sites
CAPscores_class$Class <- as.factor(CAPscores_class$Class)
CAPscores_site$Site <- as.factor(CAPscores_site$Site)


# make ordination plot for CAP model for Class
CAPplot_Class <- ggplot(CAPscores_class, aes(x=CAP1, y=CAP2)) + #sets up the plot
  geom_point(aes(CAP1, CAP2, colour = Class, shape=Class, fill= Class), size=2, stroke=1)+
  scale_shape_manual(values=c(21, 22, 23), labels=c("Background", "Occupied \nBurrow", "Unoccupied \nBurrow"))+
  scale_fill_discrete_qualitative(palette = "Dark 3", labels=c("Background", "Occupied \nBurrow", "Unoccupied \nBurrow"), alpha=0.5, nmax=6, order = c(1,3,5)) +
  scale_color_discrete_qualitative(palette = "Dark 3",labels=c("Background", "Occupied \nBurrow", "Unoccupied \nBurrow"), nmax=6, order = c(1,3,5)) +
  #scale_color_manual(values=colors30, labels=c("Background", "Occupied Burrow", "Unoccupied Burrow"))+
  #scale_fill_manual(values=alpha(colors30, 0.3), labels=c("Background", "Occupied Burrow", "Unoccupied Burrow"))+
  theme(axis.title.x = element_text(size=12, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text( size=12, family = "Arial"), 
        axis.text.y = element_blank(), axis.ticks=element_blank())+
  #theme(legend.title=element_blank(), legend.text = element_text(size=10))+
  theme(legend.position="bottom", 
        legend.title= element_blank(), 
        legend.text = element_text(size=11, margin= margin(l=1), family = "Arial"), legend.key.spacing.x = unit(2, "pt")) +
  theme(panel.border = element_rect(fill=NA, colour = "black", linewidth=1, linetype = "solid"), panel.background = element_rect(fill="white"))

CAPplot_Class

# make ordination plot for CAP model for Sites
CAPplot_Site <- ggplot(CAPscores_site, aes(x=CAP1, y=CAP2)) + #sets up the plot
  geom_point(aes(CAP1, CAP2, colour = Site, shape = Site, fill= Site), size=2, stroke=1)+
  scale_shape_manual(values=c(21, 22, 23), labels=c("Site 1", "Site 2", "Site 3"))+
  scale_fill_discrete_qualitative(palette = "Dark 3", labels=c("Site 1", "Site 2", "Site 3"), alpha=0.5, nmax=6, order = c(2,4,6)) +
  scale_color_discrete_qualitative(palette = "Dark 3", labels=c("Site 1", "Site 2", "Site 3"), nmax=6, order = c(2,4,6)) +
  #scale_color_manual(values=colors30, labels=c("Site 1", "Site 2", "Site 3"))+
  #scale_fill_manual(values=alpha(colors30, 0.3), labels=c("Site 1", "Site 2", "Site 3"))+
  theme(axis.title.x = element_text(size=12, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=12, family = "Arial"), axis.text.y = element_blank(), axis.ticks=element_blank())+
  #theme(legend.title=element_blank(), legend.text = element_text(size=10))+
  theme(legend.position="bottom", legend.title=element_blank(), 
        legend.text = element_text(size=11, margin= margin(l=1), family = "Arial"), 
        legend.key.spacing.x = unit(2, "pt"))+
  theme(panel.border = element_rect(fill=NA, colour = "black", linewidth=1, linetype = "solid"), panel.background = element_rect(fill="white"))

CAPplot_Site

# put two plots together
CAPplot_Site + CAPplot_Class + plot_annotation(tag_levels="A") & theme(plot.tag = element_text(face = 'bold', family = "Arial"))
# this is Figure 2 in the corresponding manuscript

# export plot
ggsave(filename = here("Figures", "CAP_Site&Class.png"), 
       width=8, height=4.5, units="in", dpi=300, device="png")


ggsave(filename = here("Figures", "Figure2.pdf"), 
       width=8, height=4.5, units="in", dpi=300, device=cairo_pdf)

### Examine correlations between the CAP model axes and the compounds
# this is another way to determine which compounds help differentiate sample types and colony sites
# this complements the ANOVA approach used in Step1 script

# For Class
CAPcorr_class <- read.csv(here("PRIMER/Exported Results", "Class CAP Axes Compound Correlations.csv"), header=T)

# identify compounds strongly correlated with first CAP axis:
CAP1_class <- CAPcorr_class %>% filter(CAP1<=-0.5 | CAP1>=0.5) %>% select(-CAP2)
CAP1_class

# identify compounds strongly correlated with second CAP axis:
CAP2_class <- CAPcorr_class %>% filter(CAP2<=-0.5 | CAP2>=0.5) %>% select(-CAP1)
CAP2_class

# For Site
CAPcorr_site <- read.csv(here("PRIMER/Exported Results", "Site CAP Axes Compound Correlations.csv"), header=T)

# identify compounds strongly correlated with first CAP axis:
CAP1_site <- CAPcorr_site %>% filter(CAP1<=-0.5 | CAP1>=0.5) %>% 
  select(-CAP2)
CAP1_site

# compounds strongly correlated with second CAP axis:
CAP2_site <- CAPcorr_site %>% filter(CAP2<=-0.5 | CAP2>=0.5) %>% 
  select(-CAP1)
CAP2_site


