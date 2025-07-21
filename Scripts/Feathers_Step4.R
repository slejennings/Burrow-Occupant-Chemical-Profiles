################################################################################################################
### Project: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird
### Script: Soil_Step1.R
### This script is step 4 of 6
### Script objective(s): organize and process chemical data from feather samples. Examine the repeatability across triplicate samples
################################################################################################################


######### Setup #############

# load packages
library(here)
library(tidyverse)
library(vegan)
library(perm)
library(writexl)
library(patchwork)

######### Import and combine data #########

# peak areas for feather compounds
feather_PA <- read.csv(here("Data", "featherpeakarea.csv"), header=T)
head(feather_PA)

# feather sample information
feather_sampleinfo <- read.csv(here("Data", "feathersampleinfo.csv"), header=T)
head(feather_sampleinfo)

# combine peak areas with sample information
feather <- left_join(feather_PA, feather_sampleinfo)
head(feather)

# standardize the peak area by dividing by the internal standard area and then the sample mass
StdPA_feather <- feather %>%
  mutate(Std_PA = (Area/IS_Area/Sample_Mass))

# convert to wide format
# convert data to wide format and set location and replicate info as rownames
StdPA_feather_wide <- StdPA_feather %>%
  mutate(Rep = str_sub(File, start=-3, end=-3),
         Band_Rep = paste(Band_Number, Rep, sep="_")) %>%
  dplyr::select(Band_Rep, RT, Std_PA) %>%
  pivot_wider(., names_from = RT, values_from = Std_PA) %>%
  column_to_rownames("Band_Rep")



# log(x+1) transformation
StdPA_feather_log <-log(StdPA_feather_wide[,]+1)# transform all columns

# make a matrix of pairwise distances using bray curtis dissimilarity
StdPA_feather_bray <- vegdist(StdPA_feather_log[,],"bray")

# convert the object back into a matrix
StdPA_feather_bray_matrix <- as.matrix(StdPA_feather_bray, labels=TRUE)

# turn upper triangle of matrix into NAs, including the diagonal
StdPA_feather_bray_matrix[upper.tri(StdPA_feather_bray_matrix, diag=TRUE)] <- NA

# convert data in lower triangle (not NAs) back into a long dataframe
feather_bray_df <- as.data.frame.table(StdPA_feather_bray_matrix, responseName="BrayCurtis") %>%
  filter(complete.cases(BrayCurtis)) %>%
  rename(ID1 = Var1, ID2 = Var2)

# split columns and label rows where the individual bird is the same
feather_separate <- feather_bray_df %>% 
  separate(ID1, c("Bird1","Rep1"), sep="_") %>%
  separate(ID2, c("Bird2","Rep2"), sep="_") %>%
  mutate(Match = if_else(Bird1==Bird2, "Same", "Different"))
# same indicates samples collected from the same bird, different indicates the samples came from different birds

# prepare for a permutation-based t-test
same_feather <- feather_separate %>%
  filter(Match=="Same") %>%
  dplyr::select(BrayCurtis)

diff_feather <- feather_separate %>%
  filter(Match=="Different") %>%
  select(BrayCurtis)

same_feather_vect <- same_feather[,1] # convert data into numeric vector for permTS
diff_feather_vect <- diff_feather[,1] # convert data into numeric vector for permTS

# perform monte carlo t-test with 9999 iterations
set.seed(739)
permTS(same_feather_vect, diff_feather_vect, alternative = "two.sided", method = "exact.mc", 
       control = permControl(nmc=9999, setSEED = FALSE, tsmethod = "abs"))
# highly significant: p < 0.0001


# find mean values to plot as lines on graph
mean_feather <- feather_separate %>%
  group_by(Match) %>% 
  summarise(Ave=mean(BrayCurtis))
# mean BC dissimilarity between two samples from same bird = 0.0789
# mean BC dissimilarity between two samples form different birds = 0.214

# make a density graph that includes the two lines
ggplot(feather_separate, aes(x = BrayCurtis, fill= Match, color = Match)) +
  geom_density(alpha=.3) +
  geom_vline(data= mean_feather, aes(xintercept=Ave, color=Match), linetype="dashed", linewidth =.75, show.legend=F)+
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
  scale_y_continuous(limits=c(0, 20), breaks=seq(0,20,2))


# find average across replicates
AveStdPA_feather <- StdPA_feather %>% 
  group_by(Band_Number, RT) %>% 
  summarize(Ave_StdPA=mean(Std_PA))

print(AveStdPA_feather, n=10) # look at first 10 rows

# import compound list
compounds <- read.csv(here("Data", "feathercompoundlist.csv"), header=T)
head(compounds)

# combine compound list with compound abundance
feathers_compounds <- compounds %>% 
  rename(RT = FeatherRT) %>% 
  left_join(AveStdPA_feather, .)

# add sample info about which analyses each sample should be used for
feathers_combined <- feather_sampleinfo %>%
  dplyr::select(Band_Number, Burrow, Pair, BurrowOverlap) %>%
  distinct() %>%
  left_join(feathers_compounds, .)

# export data frame for use in subsequent steps
saveRDS(feathers_combined, here("Outputs", "feathers_combined.rds"))

# confirm every bird has 155 compounds
feathers_count <- feathers_combined %>% group_by(Band_Number) %>% count()
print(feathers_count, n=56)

# how many individuals are in the data?
length(unique(feathers_count$Band_Number)) # 56

# confirm every compound has 56 samples (birds)
comp_count <- feathers_combined %>% group_by(RT) %>% count()
print(comp_count, n=155)

# get data for pairs analysis
pairs <- feathers_combined %>% filter(Pair == "Y") %>% select(-In_Burrow, -SoilRT, -BurrowOverlap)
length(unique(pairs$Burrow)) # 22 breeding pairs (44 individuals) for analysis

# export data frame for use in subsequent steps
saveRDS(pairs, here("Outputs", "pairs.rds"))

# get data for occupant and burrow analysis
occupant <- feathers_combined %>% filter(BurrowOverlap == "Y") %>% select(-Pair)
length(unique(occupant$Band_Number)) # 22 birds with soil samples from their burrow for analysis

# export data frame for use in subsequent steps
saveRDS(occupant, here("Outputs", "occupant.rds"))
