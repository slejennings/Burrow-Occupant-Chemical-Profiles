################################################################################################################
### Project: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird
### Script: Soil_Step1.R
### This script is step 5 of 6
### Script objective(s): examine the chemical profiles of birds that are part of a breeding pair 
################################################################################################################

######### Setup #############

# load packages
library(here)
library(tidyverse)
library(vegan)
library(perm)
library(writexl)
library(patchwork)

# Import pairs data
pairs <- readRDS(here("Outputs", "pairs.rds"))

############### PAIRS #######################

# subset data into different groups of compounds
# 5 groups: all compounds, plant-derived compounds, contaminants, bird-derived compounds, and bird-derived compounds that are elevated in occupied burrows
pairs_plant <- pairs %>% filter(Type=="Plant") 
pairs_contaminant <- pairs %>% filter(Type=="Contaminant") 
pairs_allbird <- pairs %>% filter(Type=="Bird") 
# elevated bird comps are: RT5.39 hexanal, RT8.23 heptanal, RT14.62 nonanal, RT17.63 decanal, RT20.47 undecanal, RT30.19 pristane
pairs_elevbird <- pairs %>% filter(RT %in% c("RT5.39", "RT8.23", "RT14.62", "RT17.63", "RT20.47", "RT30.19"))

# pivot data frames wider
pairs_wide <- pairs %>% select(-Type, -Pair, -Compound_Name) %>%
  pivot_wider(., id_cols = c(Burrow, Band_Number), names_from = RT, values_from = Ave_StdPA) %>%
  column_to_rownames(., var="Band_Number")

pairs_plant_wide <- pairs_plant %>% select(-Type, -Pair, -Compound_Name) %>%
  pivot_wider(., id_cols = c(Burrow, Band_Number), names_from = RT, values_from = Ave_StdPA) %>%
  column_to_rownames(., var="Band_Number")

pairs_contaminant_wide <- pairs_contaminant %>% select(-Type, -Pair, -Compound_Name) %>%
  pivot_wider(., id_cols = c(Burrow, Band_Number), names_from = RT, values_from = Ave_StdPA) %>%
  column_to_rownames(., var="Band_Number")

pairs_allbird_wide <- pairs_allbird %>% select(-Type, -Pair, -Compound_Name) %>%
  pivot_wider(., id_cols = c(Burrow, Band_Number), names_from = RT, values_from = Ave_StdPA) %>%
  column_to_rownames(., var="Band_Number")

pairs_elevbird_wide <- pairs_elevbird %>% select(-Type, -Pair, -Compound_Name) %>%
  pivot_wider(., id_cols = c(Burrow, Band_Number), names_from = RT, values_from = Ave_StdPA) %>%
  column_to_rownames(., var="Band_Number")

# export data frames as Excel files to be analyzed in PRIMER
#writexl::write_xlsx(pairs_wide, here("PRIMER/Imported Files", "Pairs_AllBirdComps.xlsx"))
#writexl::write_xlsx(pairs_plant_wide, here("PRIMER/Imported Files", "Pairs_PlantComps.xlsx"))
#writexl::write_xlsx(pairs_contaminant_wide, here("PRIMER/Imported Files", "Pairs_Contaminants.xlsx"))
#writexl::write_xlsx(pairs_allbird_wide, here("PRIMER/Imported Files", "Pairs_AllBirdComps.xlsx"))
#writexl::write_xlsx(pairs_elevbird_wide, here("PRIMER/Imported Files", "Pairs_ElevBirdComps.xlsx"))

# log(x+1) transform each data frame
# drop Burrow column during transformation
pairs_wide_log <- log(pairs_wide[,-1]+1) 
pairs_plant_wide_log <- log(pairs_plant_wide[,-1]+1) 
pairs_contaminant_wide_log <- log(pairs_contaminant_wide[,-1]+1) 
pairs_allbird_wide_log <- log(pairs_allbird_wide[,-1]+1) 
pairs_elevbird_wide_log <- log(pairs_elevbird_wide[,-1]+1) 

# make Bray Curtis dissimilarity matrix from each wide data frame
pairs_bray <- vegdist(pairs_wide_log,"bray")
pairs_plant_bray <- vegdist(pairs_plant_wide_log,"bray")
pairs_contaminant_bray <- vegdist(pairs_contaminant_wide_log,"bray")
pairs_allbird_bray <- vegdist(pairs_allbird_wide_log,"bray")
pairs_elevbird_bray <- vegdist(pairs_elevbird_wide_log,"bray")

################################################
# Use permanova to test for an effect of pair
# these models were also run in PRIMER to get accurate SS and MS values
# however, for one-way designs that are balanced like these models, adonis2 produces the same F-ratio and p-value as PRIMER

# get 5 random numbers to be used in set.seed() for each model
set.seed(321)
seednums <- sample(100:999, 5, replace=F)

# Burrow is the variable that groups individuals into pairs
set.seed(seednums[1])
adonis2(pairs_bray ~ as.factor(Burrow), data = pairs_wide, permutations=9999)

set.seed(seednums[2])
adonis2(pairs_plant_bray ~ as.factor(Burrow), data= pairs_plant_wide, permutations=9999)

set.seed(seednums[3])
adonis2(pairs_contaminant_bray ~ as.factor(Burrow), data= pairs_contaminant_wide, permutations=9999)

set.seed(seednums[4])
adonis2(pairs_allbird_bray ~ as.factor(Burrow), data= pairs_allbird_wide, permutations=9999)

set.seed(seednums[5])
adonis2(pairs_elevbird_bray ~ as.factor(Burrow), data= pairs_elevbird_wide, permutations=9999)


# make nMDS plots to visualize the relationship between the chemical profiles of birds
# give individuals in breeding pairs the same color and shape in these plots

# set colors for plots
# we need 22 colors for the 22 pairs
colors22<-c("#ff7690","#61c100", "#8934e4", "#008525", "#201eb2", "#ee8a00", "#0158bf", "#9cb553", "#b800ba","#27c09b", "#ff2fbd",  
            "#677100","#fa7ee1", "#6d3c00","#00aaf8","#b94600", "#263871", "#ba0033", "#df92cc", "#e09b77", "#5c255f", "#91516b")



# run nMDS for each group of compounds
# using k=3 for each here as the stress value is better

# for all compounds
pairs_nmds <- metaMDS(pairs_bray, k=3) # stress = 0.10
pairs_nmds
stressplot(pairs_nmds) # look at stress plot

# for plant compounds
pairs_plant_nmds <- metaMDS(pairs_plant_bray, k=3) # stress = 0.14
pairs_plant_nmds

# for contaminants
pairs_contaminant_nmds <- metaMDS(pairs_contaminant_bray, k=3) # stress = 0.10
pairs_contaminant_nmds

# for all bird compounds
pairs_allbird_nmds <- metaMDS(pairs_allbird_bray, k=3) # stress = 0.07
pairs_allbird_nmds

# for elevated bird compounds
pairs_elevbird_nmds <- metaMDS(pairs_elevbird_bray, k=3) # stress = 0.09
pairs_elevbird_nmds

# extract scores of individuals to data frame so they can be plotted
pairs_scores_allcomps <- as.data.frame(scores(pairs_nmds, display="sites")) %>% 
  mutate(Burrow = pairs_wide$Burrow) 

pairs_scores_plant <- as.data.frame(scores(pairs_plant_nmds, display="sites")) %>% 
  mutate(Burrow = pairs_plant_wide$Burrow) 

pairs_scores_contaminant <- as.data.frame(scores(pairs_contaminant_nmds, display="sites")) %>% 
  mutate(Burrow = pairs_contaminant_wide$Burrow) 

pairs_scores_allbird <- as.data.frame(scores(pairs_allbird_nmds, display="sites")) %>% 
  mutate(Burrow = pairs_allbird_wide$Burrow) 

pairs_scores_elevbird <- as.data.frame(scores(pairs_elevbird_nmds, display="sites")) %>% 
  mutate(Burrow = pairs_elevbird_wide$Burrow) 

# create nMDS plots

nmds.plot.allcomps <- ggplot(pairs_scores_allcomps, aes(x=NMDS1, y=NMDS2))+ # sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11),
        axis.text.y = element_blank(), 
        axis.ticks=element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1 , linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "All compounds") +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  annotate("text", x =-0.25, y =-0.2, label ="stress = 0.10" )

nmds.plot.allcomps

nmds.plot.plant <- ggplot(pairs_scores_plant, aes(x=NMDS1, y=NMDS2))+ # sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11),
        axis.text.y = element_blank(), 
        axis.ticks=element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Plant compounds") +
  theme(plot.title = element_text(hjust = 0.5, size =11)) +
  annotate("text", x =-0.6, y =-0.6, label ="stress = 0.14")

nmds.plot.plant

nmds.plot.contaminant <- ggplot(pairs_scores_contaminant, aes(x=NMDS1, y=NMDS2))+ # sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke= 1)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11),
        axis.text.y = element_blank(), 
        axis.ticks=element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Contaminant compounds") +
  theme(plot.title = element_text(hjust = 0.5, size =11)) +
  annotate("text", x =-0.2, y =-0.2, label ="stress = 0.10")

nmds.plot.contaminant

nmds.plot.allbird <- ggplot(pairs_scores_allbird, aes(x=NMDS1, y=NMDS2))+ # sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11),
        axis.text.y = element_blank(), 
        axis.ticks=element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "All bird compounds") +
  theme(plot.title = element_text(hjust = 0.5, size =11)) +
  annotate("text", x =-0.29, y =-0.32, label ="stress = 0.07")

nmds.plot.allbird

nmds.plot.elevbird <- ggplot(pairs_scores_elevbird, aes(x=NMDS1, y=NMDS2))+ # sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=12),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=12),
        axis.text.y = element_blank(), 
        axis.ticks=element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Bird compounds that were \n  elevated in occupied burrows") +
  theme(plot.title = element_text(hjust = 0.5, size =11)) +
  annotate("text", x =-0.1, y =-0.12, label ="stress = 0.09")

nmds.plot.elevbird

# combine plots. this is Figure 5 in the corresponding manuscript
nmds.plot.plant + nmds.plot.elevbird + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here("Figures", "Pairs_nMDS_PlantElevBird.png"), 
       width=8.5, height= 3.75, units="in", dpi=300, device="png")

# combine other plots for supplement
nmds.plot.allcomps + nmds.plot.contaminant + nmds.plot.allbird + 
  plot_layout(ncol=2) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here("Figures", "Pairs_nMDS_All_Bird_Contaminants.png"), 
       width=8.5, height= 6.5, units="in", dpi=300, device="png")
