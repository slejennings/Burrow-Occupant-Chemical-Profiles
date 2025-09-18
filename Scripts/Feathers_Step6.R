################################################################################################################
### Project: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird
### Script: Soil_Step1.R
### This script is step 6 of 6
### Script objective(s): compare chemical profiles of burrows (soil samples) and their avian occupants (feather samples) 
### Script objective(s) continued: examine evidence to support a transfer of chemical fingerprints from the bird to the burrow and vice versa
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
occupant <- readRDS(here("Outputs", "occupant.rds"))

# Import compound list
compounds <- read.csv(here("Data", "feathercompoundlist.csv"), header=T)

############# OCCUPANTS #################

# reduce to 75 compounds that were also detected in burrow
occupant_75 <- occupant %>% filter(In_Burrow == "Y") %>%
  rename(FeatherRT = RT)

StdPA_Soil <- readRDS(here("Outputs", "StdPA_Soil.rds"))
head(StdPA_Soil)

# put identifying info for each sample into a data frame
Soil_ID <- StdPA_Soil %>% 
  select(RT, Type, Location, Burrow, Feathers, Occupied, Pair) %>% distinct()

# find average standardized peak area for each compounds across triplicate samples
AveStdPA_Soil <- StdPA_Soil %>% 
  group_by(Location, RT) %>% 
  summarize(Ave_StdPA = mean(Std_PA)) %>%
  left_join(., Soil_ID) # combine with identifying info

# reduce to only samples from inside burrows
AveStdPA_Burrow <- AveStdPA_Soil %>% filter(Type == "Burr") %>%
  select(-Type)

# reduce to 75 compounds also detected in birds
# keep samples where we have feathers from an occupant
comps_75 <- compounds %>% filter(In_Burrow == "Y") %>% # get list of 75 compounds to keep
  rename(RT = SoilRT) %>% select(-In_Burrow)

burrow_75 <- left_join(comps_75, AveStdPA_Burrow, by="RT") %>% # use compound list to simplify burrow data
  filter(Feathers == "Y", Occupied == "Y")

# are the number of burrows and number of occupants remaining equal?
length(unique(burrow_75$Burrow)) == length(unique(occupant_75$Band_Number))

# do the burrow numbers match in both datasets (is there an occupant for each burrow?) ?
sort(unique(occupant_75$Burrow)) == sort(unique(burrow_75$Burrow))

# make burrow data wide
burrow_75_wide <- burrow_75 %>% 
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

burrow_plant_wide <- burrow_75 %>% 
  filter(Type == "Plant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

burrow_contaminant_wide <- burrow_75 %>% 
  filter(Type == "Contaminant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

burrow_allbird_wide <- burrow_75 %>% 
  filter(Type == "Bird") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

burrow_elevbird_wide <- burrow_75 %>% 
  filter(FeatherRT %in% c("RT5.39", "RT8.23", "RT14.62", "RT17.63", "RT20.47", "RT30.19")) %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# make occupant data wide
occupant_75_wide <- occupant_75 %>% 
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

occupant_plant_wide <- occupant_75 %>% 
  filter(Type == "Plant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

occupant_contaminant_wide <- occupant_75 %>% 
  filter(Type == "Contaminant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

occupant_allbird_wide <- occupant_75 %>% 
  filter(Type == "Bird") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

occupant_elevbird_wide <- occupant_75 %>% 
  filter(FeatherRT %in% c("RT5.39", "RT8.23", "RT14.62", "RT17.63", "RT20.47", "RT30.19")) %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# perform some checks to confirm the corresponding data frames for occupants and burrows have the same columns
ncol(occupant_75_wide) == ncol(burrow_75_wide)
ncol(occupant_plant_wide) == ncol(burrow_plant_wide)
ncol(occupant_contaminant_wide) == ncol(burrow_contaminant_wide)
ncol(occupant_allbird_wide) == ncol(burrow_allbird_wide)
ncol(occupant_elevbird_wide) == ncol(burrow_elevbird_wide)

# log(x+1) transform data

# burrow data
burrow_75_log <- log(burrow_75_wide+1)
burrow_plant_log <- log(burrow_plant_wide+1)
burrow_contaminant_log <- log(burrow_contaminant_wide+1)
burrow_allbird_log <- log(burrow_allbird_wide+1)
burrow_elevbird_log <- log(burrow_elevbird_wide+1)

# occupant data
occupant_75_log <- log(occupant_75_wide+1)
occupant_plant_log <- log(occupant_plant_wide+1)
occupant_contaminant_log <- log(occupant_contaminant_wide+1)
occupant_allbird_log <- log(occupant_allbird_wide+1)
occupant_elevbird_log <- log(occupant_elevbird_wide+1)


# create bray curtis dissimilarity matrices
burrow_75_bray <- vegdist(burrow_75_log, "bray")
burrow_plant_bray <- vegdist(burrow_plant_log, "bray")
burrow_contaminant_bray <- vegdist(burrow_contaminant_log, "bray")
burrow_allbird_bray <- vegdist(burrow_allbird_log, "bray")
burrow_elevbird_bray <- vegdist(burrow_elevbird_log, "bray")

occupant_75_bray <- vegdist(occupant_75_log, "bray")
occupant_plant_bray <- vegdist(occupant_plant_log, "bray")
occupant_contaminant_bray <- vegdist(occupant_contaminant_log, "bray")
occupant_allbird_bray <- vegdist(occupant_allbird_log, "bray")
occupant_elevbird_bray <- vegdist(occupant_elevbird_log, "bray")


# run procrustes analysis

# generate some random numbers for set.seed for permutation tests
set.seed(987)
seednums2 <- sample(100:999, 5, replace=F)

# all 75 compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums2[1])
(pro_test_75 <-protest(burrow_75_bray, occupant_75_bray, scores="sites", permutations=9999))

# run nMDS
burr_75_nmds <- metaMDS(burrow_75_bray, k=3)
occ_75_nmds <- metaMDS(occupant_75_bray, k=3)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_nmds_75 <- plot(procrustes(burr_75_nmds, occ_75_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_75 <- data.frame(procrustes_nmds_75$heads) %>% rownames_to_column(., var="Burrow")
points_75 <- data.frame(procrustes_nmds_75$points) %>% rownames_to_column(., var="Burrow")
plotdat_75 <- left_join(heads_75, points_75)

# set colors for plots
# we need 22 colors for the 22 pairs
colors22<-c("#ff7690","#61c100", "#8934e4", "#008525", "#201eb2", "#ee8a00", "#0158bf", "#9cb553", "#b800ba","#27c09b", "#ff2fbd",  
            "#677100","#fa7ee1", "#6d3c00","#00aaf8","#b94600", "#263871", "#ba0033", "#df92cc", "#e09b77", "#5c255f", "#91516b")


procrustes.plot_75 <- ggplot()+ 
  geom_point(data = heads_75, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_75, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_75, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, , alpha = 0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  xlim(-0.5,0.25) + ylim(-0.25, 0.25) +
  labs(title = "All compounds") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_75


# plant compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums2[2])
(pro_test_plant <-protest(burrow_plant_bray, occupant_plant_bray, scores="sites", permutations=9999))

# run nMDS
burr_plant_nmds <- metaMDS(burrow_plant_bray, k=3)
occ_plant_nmds <- metaMDS(occupant_plant_bray, k=3)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_nmds_plant <- procrustes(burr_plant_nmds, occ_plant_nmds, symmetric = TRUE, scores = "sites", scale=TRUE)
procrustes_nmds_plant <- plot(procrustes(burr_plant_nmds, occ_plant_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_plant <- data.frame(procrustes_nmds_plant$heads) %>% rownames_to_column(., var="Burrow")
points_plant <- data.frame(procrustes_nmds_plant$points) %>% rownames_to_column(., var="Burrow")
plotdat_plant <- left_join(heads_plant, points_plant)

procrustes.plot_plant <- ggplot()+ 
  geom_point(data = heads_plant, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_plant, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_plant, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, alpha = 0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Plant compounds") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_plant



# contaminants
# perform procrustes permutation-based test of significance 
set.seed(seednums2[3])
(pro_test_contaminant <-protest(burrow_contaminant_bray, occupant_contaminant_bray, scores="sites", permutations=9999))

# run nMDS
burr_contaminant_nmds <- metaMDS(burrow_contaminant_bray, k=2)
occ_contaminant_nmds <- metaMDS(occupant_contaminant_bray, k=2)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_nmds_contaminant <- procrustes(burr_contaminant_nmds, occ_contaminant_nmds, symmetric = TRUE, scores = "sites", scale=TRUE)
procrustes_nmds_contaminant <- plot(procrustes(burr_contaminant_nmds, occ_contaminant_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_contaminant <- data.frame(procrustes_nmds_contaminant$heads) %>% rownames_to_column(., var="Burrow")
points_contaminant <- data.frame(procrustes_nmds_contaminant$points) %>% rownames_to_column(., var="Burrow")
plotdat_contaminant <- left_join(heads_contaminant, points_contaminant)

procrustes.plot_contaminant <- ggplot()+ 
  geom_point(data = heads_contaminant, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_contaminant, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_contaminant, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, alpha=0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Contaminant compounds") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_contaminant

# all bird compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums2[4])
(pro_test_allbird <-protest(burrow_allbird_bray, occupant_allbird_bray, scores="sites", permutations=9999))

# run nMDS
burr_allbird_nmds <- metaMDS(burrow_allbird_bray, k=3)
occ_allbird_nmds <- metaMDS(occupant_allbird_bray, k=3)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_nmds_allbird <- procrustes(burr_allbird_nmds, occ_allbird_nmds, symmetric = TRUE, scores = "sites", scale=TRUE)
procrustes_nmds_allbird <- plot(procrustes(burr_allbird_nmds, occ_allbird_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_allbird <- data.frame(procrustes_nmds_allbird$heads) %>% rownames_to_column(., var="Burrow")
points_allbird <- data.frame(procrustes_nmds_allbird$points) %>% rownames_to_column(., var="Burrow")
plotdat_allbird <- left_join(heads_allbird, points_allbird)

procrustes.plot_allbird <- ggplot()+ 
  geom_point(data = heads_allbird, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_allbird, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_allbird, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, alpha = 0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "All bird compounds") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_allbird

# elevated bird compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums2[5])
(pro_test_elevbird <-protest(burrow_elevbird_bray, occupant_elevbird_bray, scores="sites", permutations=9999))

# run nMDS
burr_elevbird_nmds <- metaMDS(burrow_elevbird_bray, k=3)
occ_elevbird_nmds <- metaMDS(occupant_elevbird_bray, k=3)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_nmds_elevbird <- procrustes(burr_elevbird_nmds, occ_elevbird_nmds, symmetric = TRUE, scores = "sites", scale=TRUE)
procrustes_nmds_elevbird <- plot(procrustes(burr_elevbird_nmds, occ_elevbird_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_elevbird <- data.frame(procrustes_nmds_elevbird$heads) %>% rownames_to_column(., var="Burrow")
points_elevbird <- data.frame(procrustes_nmds_elevbird$points) %>% rownames_to_column(., var="Burrow")
plotdat_elevbird <- left_join(heads_elevbird, points_elevbird)

procrustes.plot_elevbird <- ggplot()+ 
  geom_point(data = heads_elevbird, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_elevbird, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_elevbird, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, alpha=0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors22)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Bird compounds that were \n   elevated in occupied burrows") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_elevbird

# combine plots for plant and elevated bird compounds
# this is figure 4 in the corresponding manuscript
procrustes.plot_plant + procrustes.plot_elevbird + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', family = "Arial"))

ggsave(filename = here("Figures", "Procrustes_PlantElevBird.png"), 
       width=8, height= 4, units="in", dpi=300, device="png")

ggsave(filename = here("Figures", "Figure4.pdf"), 
       width=8, height= 4, units="in", dpi=300, device=cairo_pdf)

# combine the remaining plots for the supplement
procrustes.plot_75 + procrustes.plot_contaminant + procrustes.plot_allbird + 
  plot_layout(ncol=2) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', family = "Arial"))

ggsave(filename = here("Figures", "Procrustes_All_Contaminant_AllBird.png"), 
       width=8, height= 6.75, units="in", dpi=300, device="png")

##################################################################################################################
##################################################################################################################
##################################################################################################################

# Here, we complete the same analyses as above but using the 10 feather samples collected 4-7 days after the burrow soil
# these individuals represent some of the mates of the birds analyzed in the previous step
# the objective is to see if the same findings hold when using birds collected at a later time point when we may expect to see a lower degree of overlap 


# Import data
feathers_combined <- readRDS(here("Outputs", "feathers_combined.rds")) # feather data
StdPA_Soil <- readRDS(here("Outputs", "StdPA_Soil.rds")) # soil data


# get the of burrows that have feather samples that came from pairs
# from this, keep the Burrow ID
burrow_pairs <- feathers_combined %>% 
  filter(Pair == "Y") %>%
  filter(BurrowOverlap == "Y") %>%
  pull(Burrow) %>% 
  unique()

length(unique(burrow_pairs)) # these are the 10 burrows we want for this analysis

# identify the birds with feather samples were not present in the burrow on the same day the soil was collected
# feather samples from these individuals were collected on a subsequent visit that occurred 4-7 days after the soil sampling
occupant_2nd <- feathers_combined %>%
  filter(Burrow %in% burrow_pairs) %>%
  filter(BurrowOverlap == "N")

# reduce the chemical profiles of these birds to the 75 compounds that were also detected in burrow
occupant_2nd_75 <- occupant_2nd %>% filter(In_Burrow == "Y") %>%
  rename(FeatherRT = RT)

# next, work to prepare the soil data
# put the identifying info for each soil sample into a data frame
Soil_ID <- StdPA_Soil %>% 
  select(RT, Type, Location, Burrow, Feathers, Occupied, Pair) %>% distinct()

# find average standardized peak area for each compounds across triplicate samples
AveStdPA_Soil <- StdPA_Soil %>% 
  group_by(Location, RT) %>% 
  summarize(Ave_StdPA = mean(Std_PA)) %>%
  left_join(., Soil_ID) # combine with identifying info

# reduce to only samples from inside burrows
AveStdPA_Burrow <- AveStdPA_Soil %>% filter(Type == "Burr") %>%
  select(-Type)

# reduce to 75 compounds also detected in birds
# keep the samples where we have feathers from an occupant
comps_75 <- compounds %>% filter(In_Burrow == "Y") %>% # get list of 75 compounds to keep
  rename(RT = SoilRT) %>% select(-In_Burrow)

# reduce to the 10 burrows with 2nd feather sample
burrow_2nd_75 <- left_join(comps_75, AveStdPA_Burrow, by="RT") %>% # use compound list to simplify burrow data
  filter(Feathers == "Y", Occupied == "Y") %>%
  filter(Burrow %in% burrow_pairs) 

# are the number of burrows and number of occupants remaining equal?
length(unique(burrow_2nd_75$Burrow)) == length(unique(occupant_2nd_75$Band_Number))

# do the burrow numbers match in both datasets (is there an occupant for each burrow?) ?
sort(unique(occupant_2nd_75$Burrow)) == sort(unique(burrow_2nd_75$Burrow))

# convert the burrow data to wide format

# all 75 compounds
burrow_2nd_75_wide <- burrow_2nd_75 %>% 
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# plant compounds
burrow_2nd_plant_wide <- burrow_2nd_75 %>% 
  filter(Type == "Plant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# contaminants
burrow_2nd_contaminant_wide <- burrow_2nd_75 %>% 
  filter(Type == "Contaminant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# all bird compounds
burrow_2nd_allbird_wide <- burrow_2nd_75 %>% 
  filter(Type == "Bird") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# six elevated bird compounds
burrow_2nd_elevbird_wide <- burrow_2nd_75 %>% 
  filter(FeatherRT %in% c("RT5.39", "RT8.23", "RT14.62", "RT17.63", "RT20.47", "RT30.19")) %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# convert the the occupant data to wide format 
# all 75 compounds
occupant_2nd_75_wide <- occupant_2nd_75 %>% 
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# plant compounds
occupant_2nd_plant_wide <- occupant_2nd_75 %>% 
  filter(Type == "Plant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# contaminants
occupant_2nd_contaminant_wide <- occupant_2nd_75 %>% 
  filter(Type == "Contaminant") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# all bird compounds
occupant_2nd_allbird_wide <- occupant_2nd_75 %>% 
  filter(Type == "Bird") %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# six elevated bird compounds
occupant_2nd_elevbird_wide <- occupant_2nd_75 %>% 
  filter(FeatherRT %in% c("RT5.39", "RT8.23", "RT14.62", "RT17.63", "RT20.47", "RT30.19")) %>%
  pivot_wider(., names_from = FeatherRT, values_from = Ave_StdPA, id_cols = Burrow) %>%
  arrange(Burrow) %>%
  column_to_rownames(., var="Burrow")

# perform some checks to confirm the corresponding data frames for occupants and burrows have the same columns
ncol(occupant_2nd_75_wide) == ncol(burrow_2nd_75_wide)
ncol(occupant_2nd_plant_wide) == ncol(burrow_2nd_plant_wide)
ncol(occupant_2nd_contaminant_wide) == ncol(burrow_2nd_contaminant_wide)
ncol(occupant_2nd_allbird_wide) == ncol(burrow_2nd_allbird_wide)
ncol(occupant_2nd_elevbird_wide) == ncol(burrow_2nd_elevbird_wide)

# log(x+1) transform data

# burrow data
burrow_2nd_75_log <- log(burrow_2nd_75_wide+1)
burrow_2nd_plant_log <- log(burrow_2nd_plant_wide+1)
burrow_2nd_contaminant_log <- log(burrow_2nd_contaminant_wide+1)
burrow_2nd_allbird_log <- log(burrow_2nd_allbird_wide+1)
burrow_2nd_elevbird_log <- log(burrow_2nd_elevbird_wide+1)

# occupant data
occupant_2nd_75_log <- log(occupant_2nd_75_wide+1)
occupant_2nd_plant_log <- log(occupant_2nd_plant_wide+1)
occupant_2nd_contaminant_log <- log(occupant_2nd_contaminant_wide+1)
occupant_2nd_allbird_log <- log(occupant_2nd_allbird_wide+1)
occupant_2nd_elevbird_log <- log(occupant_2nd_elevbird_wide+1)

# create bray curtis dissimilarity matrices
# burrow data
burrow_2nd_75_bray <- vegdist(burrow_2nd_75_log, "bray")
burrow_2nd_plant_bray <- vegdist(burrow_2nd_plant_log, "bray")
burrow_2nd_contaminant_bray <- vegdist(burrow_2nd_contaminant_log, "bray")
burrow_2nd_allbird_bray <- vegdist(burrow_2nd_allbird_log, "bray")
burrow_2nd_elevbird_bray <- vegdist(burrow_2nd_elevbird_log, "bray")

# occupant data
occupant_2nd_75_bray <- vegdist(occupant_2nd_75_log, "bray")
occupant_2nd_plant_bray <- vegdist(occupant_2nd_plant_log, "bray")
occupant_2nd_contaminant_bray <- vegdist(occupant_2nd_contaminant_log, "bray")
occupant_2nd_allbird_bray <- vegdist(occupant_2nd_allbird_log, "bray")
occupant_2nd_elevbird_bray <- vegdist(occupant_2nd_elevbird_log, "bray")


# run Procrustes analysis

# generate some random numbers for set.seed for permutation tests
set.seed(425)
seednums3 <- sample(100:999, 5, replace=F)

# all 75 compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums3[1])
(pro_test_2nd_75 <-protest(burrow_2nd_75_bray, occupant_2nd_75_bray, scores="sites", permutations=9999))

# contaminants
# perform procrustes permutation-based test of significance 
set.seed(seednums3[3])
(pro_test_2nd_contaminant <-protest(burrow_2nd_contaminant_bray, occupant_2nd_contaminant_bray, scores="sites", permutations=9999))

# all bird compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums3[4])
(pro_test_2nd_allbird <-protest(burrow_2nd_allbird_bray, occupant_2nd_allbird_bray, scores="sites", permutations=9999))

# plant compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums3[2])
(pro_test_2nd_plant <-protest(burrow_2nd_plant_bray, occupant_2nd_plant_bray, scores="sites", permutations=9999))

# this finding is consistent with the original analysis - make a plot to visualize the relationships
# run nMDS - due to small sample sizes, warning messages may appear
burr_2nd_plant_nmds <- metaMDS(burrow_2nd_plant_bray, k=3)
occ_2nd_plant_nmds <- metaMDS(occupant_2nd_plant_bray, k=3)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_2nd_nmds_plant <- procrustes(burr_2nd_plant_nmds, occ_2nd_plant_nmds, symmetric = TRUE, scores = "sites", scale=TRUE)
procrustes_2nd_nmds_plant <- plot(procrustes(burr_2nd_plant_nmds, occ_2nd_plant_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_2nd_plant <- data.frame(procrustes_2nd_nmds_plant$heads) %>% rownames_to_column(., var="Burrow")
points_2nd_plant <- data.frame(procrustes_2nd_nmds_plant$points) %>% rownames_to_column(., var="Burrow")
plotdat_2nd_plant <- left_join(heads_2nd_plant, points_2nd_plant)

# set colors for plots
# we need 10 colors for the 10 birds/burrows
colors10<-c("#ff7690","#61c100", "#8934e4", "#008525", "#201eb2", "#ee8a00", "#0158bf", "#9cb553", "#b800ba","#27c09b")

procrustes.plot_2nd_plant <- ggplot()+ 
  geom_point(data = heads_2nd_plant, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_2nd_plant, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_2nd_plant, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, alpha = 0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors10)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Plant compounds") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_2nd_plant


# elevated bird compounds
# perform procrustes permutation-based test of significance 
set.seed(seednums3[5])
(pro_test_2nd_elevbird <-protest(burrow_2nd_elevbird_bray, occupant_2nd_elevbird_bray, scores="sites", permutations=9999))

# this finding is consistent with the original analysis - make a plot to visualize the relationships
# run nMDS - due to small sample sizes, warning messages may appear
burr_2nd_elevbird_nmds <- metaMDS(burrow_2nd_elevbird_bray, k=3)
occ_2nd_elevbird_nmds <- metaMDS(occupant_2nd_elevbird_bray, k=3)

# use procrustes() to rotate and scale two nMDS plots for plotting
procrustes_nmds_2nd_elevbird <- procrustes(burr_2nd_elevbird_nmds, occ_2nd_elevbird_nmds, symmetric = TRUE, scores = "sites", scale=TRUE)
procrustes_nmds_2nd_elevbird <- plot(procrustes(burr_2nd_elevbird_nmds, occ_2nd_elevbird_nmds, symmetric = TRUE, scores = "sites", scale=TRUE))

# pull coordinates for plotting
heads_2nd_elevbird <- data.frame(procrustes_nmds_2nd_elevbird$heads) %>% rownames_to_column(., var="Burrow")
points_2nd_elevbird <- data.frame(procrustes_nmds_2nd_elevbird$points) %>% rownames_to_column(., var="Burrow")
plotdat_2nd_elevbird <- left_join(heads_2nd_elevbird, points_2nd_elevbird)

procrustes.plot_2nd_elevbird <- ggplot()+ 
  geom_point(data = heads_2nd_elevbird, aes(NMDS1, NMDS2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_point(data = points_2nd_elevbird, aes(X1, X2, colour = factor(Burrow), shape = factor(Burrow)), size = 2, stroke=1)+
  geom_segment(data = plotdat_2nd_elevbird, aes(x=NMDS1, y=NMDS2, xend=X1, yend=X2, colour = factor(Burrow)), linewidth = 1, alpha=0.6)+
  scale_shape_manual(values = rep(c(0,1,2,3,4,5,6,7,8,9,10), times=2))+
  scale_color_manual(values = colors10)+
  theme(axis.title.x = element_text(size=11, family = "Arial"),
        axis.text.x  = element_blank(), 
        axis.title.y = element_text(size=11, family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="none") +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=1, linetype = "solid"), 
        panel.background = element_rect(fill = "white")) +
  labs(title = "Bird compounds that were \n   elevated in occupied burrows") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, family = "Arial"))
procrustes.plot_2nd_elevbird

# combine plots for plant and elevated bird compounds into a single panel for supplementary materials
procrustes.plot_2nd_plant + procrustes.plot_2nd_elevbird + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', family = "Arial"))

ggsave(filename = here("Figures", "Procrustes_PlantElevBird_Supplementary.png"), 
       width=8, height= 4, units="in", dpi=300, device="png")
