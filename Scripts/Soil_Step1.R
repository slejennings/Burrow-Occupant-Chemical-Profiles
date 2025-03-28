################################################################################################################
### Project: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird
### Script: Soil_Step1.R
### This script is step 1 of 6
### Script objective(s): process soil data, example repeatability of triplicate samples (subsamples)
################################################################################################################

######### Setup #############

# load packages
library(here)
library(tidyverse)
library(broom)
library(vegan)
library(perm)
library(RColorBrewer)
library(ecodist)
library(writexl)
library(geodist) 
library(patchwork)

######### Import and combine data #########

# peak areas for soil compounds
soil_PA <- read.csv(here("Data", "soilpeakarea.csv"), header=T)
head(soil_PA)

# soil sample information
soil_sampleinfo <- read.csv(here("Data", "soilsampleinfo.csv"), header=T)
head(soil_sampleinfo)

# combine peak areas with sample information
soil <- left_join(soil_PA, soil_sampleinfo)
head(soil)

# get the dry soil mass
# the value for water percentage is the percent of the starting mass that remained after drying
# therefore, larger values are from drier burrows. 
# if the value is 0.75, this means that 25% of the weighed soil was water and 75% of it was organic (non-water) material
# water percentage to get adjusted mass
soil_dry <- soil %>% mutate(Adjust_Mass=(Wet_Soil_Mass*WaterPercentage))

# standardize the peak area by dividing by the internal standard area and then the adjusted (dry) soil mass
StdPA_Soil <- soil_dry %>%
  mutate(Std_PA = (Area/IS_Area/Adjust_Mass))

head(StdPA_Soil) 

# export this data frame for use in subsequent steps
saveRDS(StdPA_Soil, here("Outputs", "StdPA_Soil.rds"))

# how many compounds per sample?
numcomps <- StdPA_Soil %>% group_by(File) %>% filter(Area >0) %>% count()
length(unique(StdPA_Soil$RT)) # number of unique compounds in soil samples
mean(numcomps$n) # average number of compounds per sample
sd(numcomps$n) # standard deviation
boxplot(numcomps$n)

###########################################################################
#### Examine repeatability of triplicate samples ####

# each background or burrow location was analyzed in triplicate, so there are 3 subsamples associated with it 
# first, we need to test whether the replicate subsamples are consistent or repeatable in their chemical profiles
# it is important to confirm this to ensure the methods are measuring real difference in the chemical profiles of different colony locations
# one way to answer this question is to compare whether the pairwise multivariate distance between two subsamples from the same location
# is lower than the pairwise multivariate distance between two subsamples from different locations
# if this is true, we have good support that the triplicate subsamples from the same location are consistent/repeatible in their chemical profiles

# convert data to wide format and set location and replicate info as row names
StdPA_Soil_wide <- StdPA_Soil %>%
  mutate(Rep = str_sub(File, start=-3, end=-3),
         Burr_Rep = paste(Location, Rep, sep="_")) %>%
  dplyr::select(Burr_Rep, RT, Std_PA) %>%
  pivot_wider(., names_from = RT, values_from = Std_PA) %>%
  column_to_rownames("Burr_Rep")

# log (x+1) transform peak areas
StdPA_Soil_log<-log(StdPA_Soil_wide[,]+1) # transform all columns

# make a matrix of pairwise distances using bray curtis dissimilarity
StdPA_Soil_bray <- vegdist(StdPA_Soil_log[,],"bray")

# convert the object back into a matrix
StdPA_Soil_bray_matrix <- as.matrix(StdPA_Soil_bray, labels=TRUE)

# turn upper triangle of matrix into NAs, including the diagonal
StdPA_Soil_bray_matrix[upper.tri(StdPA_Soil_bray_matrix, diag=TRUE)] <- NA

# convert data in lower triangle (not NAs) back into a long data frame
soil_bray_df <- as.data.frame.table(StdPA_Soil_bray_matrix, responseName="BrayCurtis") %>%
  filter(complete.cases(BrayCurtis)) %>%
  rename(ID1 = Var1, ID2 = Var2)

# split columns and label rows where the location is the same
soil_separate <- soil_bray_df %>% 
  separate(ID1, c("Loc1","Rep1"), sep="_") %>%
  separate(ID2, c("Loc2","Rep2"), sep="_") %>%
  mutate(Match = if_else(Loc1==Loc2, "Same", "Different"))
# same indicates subsamples collected from the same location, different indicates the subsamples came from different locations

# prepare for a permutation-based t-test
same_soil <- soil_separate %>%
  filter(Match=="Same") %>%
  dplyr::select(BrayCurtis)

diff_soil <- soil_separate %>%
  filter(Match=="Different") %>%
  select(BrayCurtis)

same_soil_vect <- same_soil[,1] # convert data into numeric vector for permTS
diff_soil_vect <- diff_soil[,1] # convert data into numeric vector for permTS

# perform monte carlo t-test with 9999 iterations
set.seed(456)
perm_ttest <- permTS(same_soil_vect, diff_soil_vect, alternative = "two.sided", method = "exact.mc", 
       control = permControl(nmc=9999, setSEED = FALSE, tsmethod = "abs"))
perm_ttest
# the pairwise bray curtis dissimilarity between two subsamples from different locations is on average 0.3 greater than the pairwise dissimilarity between two subsamples from the same location
# bray curtis takes on value of 0 to 1
# highly significant: p<0.001

# save the model
saveRDS(perm_ttest, here("Models", "perm_ttest.rds"))

# find mean values to plot as lines on graph
mean_soil <- soil_separate %>%
  group_by(Match) %>% 
  summarise(Ave=mean(BrayCurtis))
# mean BC dissimilarity between two samples from same location = 0.108
# mean BC dissimilarity between two samples form different locations = 0.404

# make a density graph that includes the two lines to visualize differences
ggplot(soil_separate, aes(x=BrayCurtis, fill=Match, color=Match)) +
  geom_density(alpha=.3) +
  geom_vline(data= mean_soil, aes(xintercept=Ave, color=Match), linetype="dashed", linewidth =.75, show.legend=F)+
  labs(x="Pairwise Bray Curtis Dissimilarity", y="Fraction of Differences")+
  theme(axis.title.x = element_text(size=12),
        axis.text.x  = element_text( size=10), 
        axis.title.y = element_text( size=12),
        axis.text.y  = element_text( size=10))+
  theme(panel.background = element_rect(fill="White"))+
  theme(axis.line = element_line( linewidth =.5, linetype = "solid"))+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  theme(legend.justification=c(1,1), 
        legend.position.inside =c(1,1), 
        legend.title=element_blank(), 
        legend.text=element_text(size=10))+
  scale_x_continuous(limits=c(0, 1))+
  scale_y_continuous(limits=c(0, 11), breaks=seq(0,10,2))
