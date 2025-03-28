################################################################################################################
### Project: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird
### Script: Soil_Step1.R
### This script is step 2 of 6
### Script objective(s): Perform qualitative and quantitative comparisons of different sample types (background vs burrow)
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

######### import data #########

# import standardized peak areas for soil compounds
StdPA_Soil <- readRDS(here("Outputs", "StdPA_Soil.rds"))
head(StdPA_Soil)

# find average and standard deviation of peak area and count of occurrence for each compound in background samples
Bkgd_Stats <-StdPA_Soil %>% 
  filter(Type == "Bkgd") %>% 
  group_by(RT) %>%
  summarize(Bkgd_Ave = mean(Std_PA), Bkgd_SD = sd(Std_PA), Bkgd_Count=sum(Std_PA>0))

# find average and standard deviation of peak area and count of occurrence for each compound in burrow samples
Burr_Stats <-StdPA_Soil %>% 
  filter(Type == "Burr") %>% 
  group_by(RT) %>%
  summarize(Burr_Ave = mean(Std_PA), Burr_SD = sd(Std_PA), Burr_Count=sum(Std_PA>0))

# combine burrow and background stats into a single data frame
StdPA_Stats <- inner_join(Bkgd_Stats, Burr_Stats, by="RT")
head(StdPA_Stats)

################################################################################################################
##### Step 1: Explore qualitative differences in the chemical composition of burrow and background samples #####

# Are any compounds entirely missing from one sample type (either burrow or background)?
Missing <- StdPA_Stats %>% filter(Burr_Count==0|Bkgd_Count==0)
nrow(Missing) # no

# Are any compounds rare in one sample type (either burrow or background) but common in the other?
# And are there any compounds that are rare in both sample types?
# there are 132 files of each type of sample (e.g., 44 background samples with 3 replicates each)
# defining rare as compounds that were detected in 25% or less of samples (n= 33)
# defining common as compounds that were detected in 75% or more of samples (n = 99)
BothRare <- StdPA_Stats %>% filter(Burr_Count<=33 & Bkgd_Count<=33) # rare in both
nrow(BothRare) # no compounds meet this criteria
BurrRare <- StdPA_Stats %>% filter(Burr_Count<=33 & Bkgd_Count>=99) # rare in burrow, common in bkgd
nrow(BurrRare) # no compounds meet this criteria
BkgdRare <- StdPA_Stats %>% filter(Burr_Count>=99 & Bkgd_Count<=33) # common in burrow, rare in bkgd
nrow(BkgdRare) # no compounds meet this criteria

# Are there any compounds that are not common (rare) in both sample types?
# Meaning, are there situations where fewer than 75% of samples in both the bkgd and burr samples were missing the compound?
BothLess99 <- StdPA_Stats %>% filter(Burr_Count<99 & Bkgd_Count<99)
nrow(BothLess99) # 23 compounds meet this criteria
print(BothLess99, n=Inf)

################################################################################################################
#### Step 2: Identify quantitative differences in the chemical composition of burrow and background samples ####

# Goal: use ANOVA to compare the mean abundance of compounds between bkgd, occupied, and unoccupied burrows
# This is to give a general sense of which compounds differ between these sample types
# multivariate analyses will also be used later as a complementary approach
# removing all the rare RTs (compounds) in BothLess99 as these seem less useful from a sensory standpoint to differentiate landscape features like burrows from background soil

RTdrop <- BothLess99$RT # make list of RTs to remove

StdPA_Class <- StdPA_Soil %>%
  mutate(Class = ifelse(Type == "Bkgd", "Bkgd",
                        ifelse(Type== "Burr" & Occupied == "Y", "Occ_Burrow", "Unocc_Burrow"))) %>% # mark occupied and unoccupied burrows
  filter(!RT %in% RTdrop) %>% # remove RTs contained in RTdrop
  dplyr::select(RT, Std_PA, Class)

head(StdPA_Class)
str(StdPA_Class)
# make categorical variables into factors
StdPA_Class$RT <- as.factor(StdPA_Class$RT)
StdPA_Class$Class <- as.factor(StdPA_Class$Class)

# perform ANOVA for each compound
anova_StdPA <- StdPA_Class %>% 
  arrange(RT) %>% 
  group_nest(RT) %>% 
  mutate(aov = map(data, ~aov(Std_PA ~ Class, data = .)),
         tidyaov = map(aov, broom::tidy)) %>%
  unnest(tidyaov) %>%
  dplyr::select(-data, -aov) %>%
  filter(term == "Class")
View(anova_StdPA)

# apply benjamini-hochberg procedure to decrease false discovery rate
anova_StdPA$p.adjust <- p.adjust(anova_StdPA$p.value, method="BH")

anova_sigP <- anova_StdPA %>%
  filter(p.adjust < 0.05) %>%
  mutate(p.adjust=round(p.adjust, 4))

View(anova_sigP)
nrow(anova_sigP) # 76 compounds

RTs_sigP <- anova_sigP$RT # make a vector containing these compounds

# perform post-hoc comparisons for models that were significant
# p-values from Tukey HSD test are already adjusted for multiple comparisons
tukey_StdPA <- StdPA_Class %>%
  arrange(RT) %>%
  filter(RT %in% RTs_sigP) %>% # reduce to compounds identified above
  group_nest (RT) %>%
  mutate(tukey = map(data, ~ broom::tidy(TukeyHSD(aov(Std_PA ~ Class, data = .))))) %>%
  unnest(tukey) %>%
  dplyr::select(-data) %>%
  mutate(sigP = ifelse(adj.p.value < 0.05, "Y", "N"))

View(tukey_StdPA)

# Three different combinations of contrasts are possible in the post-hoc comparisons:
# Occ_Burrow-Bkgd
# Unocc_Burrow-Bkgd
# Unocc_Burrow-Occ_Burrow
# each will either be labeled as Y for significant difference, or N for no difference between two sample types

# look at the combinations of contrasts that exist in the data
contrast <-
  tukey_StdPA %>%
  dplyr::select(RT, contrast, sigP) %>%
  pivot_wider(., id_cols=RT, names_from=contrast, values_from=sigP) %>%
  rowwise() %>% mutate(ContrastResult = paste(`Occ_Burrow-Bkgd`, `Unocc_Burrow-Bkgd`, `Unocc_Burrow-Occ_Burrow`, sep=" ")) %>%
  dplyr::select(RT, ContrastResult)

# get estimated differences in average peak area
contrast_est <- tukey_StdPA %>%
  dplyr::select(RT, contrast, estimate) %>%
  pivot_wider(., id_cols=RT, names_from=contrast, values_from=estimate) %>%
  inner_join(., contrast)

print(contrast_est, n=Inf)       

# what are the unique combinations observed for ContrastResult?
unique(contrast_est$ContrastResult)
# "YNY" "YYN" "YNN" "NYY" "NYN" "NNY"

contrast_est %>% count(ContrastResult) # look at number of each

# Y Y N = indicates a compound that differs between the background and burrows of either kind (occ and unocc burrows) 
# these compounds could possibly be related to an above vs underground difference
YYN <- contrast_est %>%
  filter(ContrastResult=="Y Y N")
print(YYN, n=Inf)
# there are 28 compounds

# Y N Y = indicates a compounds that differentiate occupied burrows from both background and unoccupied burrows
# this result is possibly related to being actively inhabited by birds and is of particular interest in relation to olfactory burrow recognition
YNY <- contrast_est %>%
  filter(ContrastResult=="Y N Y")
YNY
# there are 10 compounds

# N Y Y = this indicates a compound in unoccupied burrows that differs from bkgd and occupied burrows (but bkgd and occupied burrows do not differ)
# this result is possibly due to less frequent disturbance to soil in unoccupied burrows than would occur in bkgd or occupied burrows
NYY <- contrast_est %>%
  filter(ContrastResult=="N Y Y")
print(NYY, n=Inf)
# there are 8 compounds

# Y N N = indicates a compound that differs between occupied burrows and the background, but no other differences 
# it is harder to come up with a clear potential biological explanation for this result
YNN <- contrast_est %>%
  filter(ContrastResult=="Y N N")
print(YNN, n=Inf)
# there are 26 compounds
# many are more elevated in background compared to occupied burrows (denoted by negative values in first column)

# N N Y and N Y N. These outcomes are much less common. Only 2 compounds each
# NNY indicates a compound that differs between occupied and unoccupied burrows but not for other comparisons
 # NYN differ for unoccupied burrows and background, but not for other comparisons
# again, harder to come up with potential biological explanations for these results
# however, these compounds will likely inform multivariate analyses

NNY <- contrast_est %>%
  filter(ContrastResult=="N N Y")
print(NNY, n=Inf)

NYN <- contrast_est %>%
  filter(ContrastResult=="N Y N")
print(NYN, n=Inf)

# create table with average peak area (abundance) for each RT and sample type
aveabundance <- StdPA_Class %>%
  arrange(RT) %>%
  filter(RT %in% RTs_sigP) %>%
  group_by(RT, Class) %>%
  summarize(AveAbund = round(mean(Std_PA),3)) %>%
  pivot_wider(., names_from=Class, values_from = AveAbund)

# the results of these steps are summarized in Table S4 in the corresponding manuscript

################################################################################################################
########### Step 3: Examine the compounds that differentiate occupied burrows from other sample types ##########
### return to YNY
# make plots for YNY as these chemicals are particular interesting because they may inform olfactory burrow recognition
occburrcomps <- YNY$RT

# prep data for plots
plotdat <- StdPA_Class %>%
  filter(RT %in% occburrcomps)

# make plots
occburr_plots <- plotdat %>%
  mutate(RTlabel = RT) %>%
  group_by(RT) %>%
  group_map(.f = ~ggplot(data=.x, aes(x = Class, y = Std_PA, fill = Class)) +
              geom_boxplot() +
              xlab("Sample Type") + 
              ylab("Standardized Peak Area") + 
              ggtitle(paste0(" ", .x$RTlabel)) + 
              theme(legend.title = element_blank(),
                    legend.position = "top",
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "darkgrey"))
  )

occburr_plots
# plots show that 8 out of 10 compounds are elevated in occupied burrows compared with other sample types
# what are these chemicals?
# 5 aldehydes. all common and abundant in feathers: 
# RT5.39 = hexanal, RT8.23 = heptanal, RT14.61 = nonanal , RT17.56 = decanal, RT20.41 = undecanal
# RT30.11 = pristane. This is the most abundant compound in feathers
# these 6 compounds are key components of the feather chemical profiles!
# The fact that they are elevated in occupied burrows suggests they are derived from the bird occupants

# 2 compounds were elevated in occupied burrows but less likely to be derived from birds
# RT8.10 = nonane. Alkane. This is not detected in feathers. But could still be related to burrow being occupied. Has multiple potential sources
# RT33.75 = homosalate. Also detected in feathers. This is a component of cosmetics (sunscreen), so likely a contaminant

# 2 compounds are lower in occupied burrows compared with other sample types
# RT10.77 is a monoterpene with m/z 107 as dominant peak in mass spectrum. Possible ID is 2-Methylenebornane. Not detected on feathers
# RT18.59 is an unidentified phthalate, so also likely a contaminant. Not detected on feathers

##########################################

# In the next step, we focus on the 6 compounds that are ubiquitous and abundant in occupied burrows and feathers
# make nMDS plot for occupied burrows replicate samples to see if they are visually differentiated using 
# a) 6 bird compounds that are elevated in occupied burrows
# b) using all compounds measured in their chemical profiles

# for the compounds that were elevated in occupied burrows, how much more abundant are they than in background or unoccupied burrow samples?
compare <- plotdat %>% filter(RT %in% c("RT5.39", "RT8.23", "RT14.61", "RT17.56", "RT20.41", "RT30.11")) %>%
  mutate(Type=if_else(Class=="Occ_Burrow", "Occ_Burrow", "Other")) %>%
  group_by(Type) %>% summarize(ave=mean(Std_PA))

compare[1,2]/compare[2,2] # divide ave abundance in occupied burows by ave abundance in bkgd and unoccupied burrows

# start with long data frame, use filter to get RTs, then pivot wider
occburr_elevbird <- StdPA_Soil %>%
  mutate(Class = ifelse(Type == "Bkgd", "Bkgd",
                        ifelse(Type== "Burr" & Occupied == "Y", "Occ_Burrow", "Unocc_Burrow"))) %>% # mark occupied and unoccupied burrows
  filter(RT %in% c("RT5.39", "RT8.23", "RT14.61", "RT17.56", "RT20.41", "RT30.11") & Class=="Occ_Burrow") %>%
  dplyr::select(File, RT, Location, Std_PA) %>%
  pivot_wider(., names_from = RT, values_from = Std_PA)

# separate id variables
occburr_elevbird_id <-occburr_elevbird %>% dplyr::select(File, Location)

# separate data, make file the rownames
colnames(occburr_elevbird)
occburr_elevbird <- column_to_rownames(occburr_elevbird, "File")

# log (x+1) transform
# drop first column which is Location when doing log transformation
occburr_elevbird_log <- log(occburr_elevbird[,-1]+1)

# use transformed data to make bray curtis distance matrix
occburr_elevbird_bray <- vegdist(occburr_elevbird_log,"bray")

# make nMDS plot
length(unique(occburr_elevbird$Location)) # how many occupied burrows are in the data? 30

# set colors for plot
# there are 30 burrows, so we need 30 unique colors
colors30 <- c( "#3e7966","#9645e0","#4dc041","#5028a3","#86b740","#d246c5","#50b97e","#585fd2","#a6a432","#b474d3","#457830","#d14b90","#5ebab0","#df462a","#6987d3","#d7a038",
               "#71327b","#95a669","#d6445c","#5b9fc5","#ce7133","#43497a","#cca270","#6c2841","#3d4a1f","#bf93c5","#845f30","#a55e75","#873224","#dc8c88")

# run nMDS 
occburr_elevbird_nmds <- metaMDS(occburr_elevbird_bray, k=2)
occburr_elevbird_nmds # stress = 0.123
stressplot(occburr_elevbird_nmds) # look at stress plot

# extract scores from nMDS and put them into a data frame for plotting
scores_occburr_elevbird <- as.data.frame(scores(occburr_elevbird_nmds, display="sites")) %>% 
  mutate(Burrow=occburr_elevbird$Location)

# create nMDS plots
NMDSplot_occburr_elevbird <- ggplot(scores_occburr_elevbird, aes(x=NMDS1, y=NMDS2)) + #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Burrow), shape=factor(Burrow)), size=2.5, stroke=1)+
  scale_shape_manual(values=rep(c(0,1,2,3,4,5,6,7,8,9,10),times=3))+
  scale_color_manual(values=colors30)+
  theme(axis.title.x = element_text(size=12),axis.text.x  = element_blank(), axis.title.y = element_text( size=12), axis.text.y = element_blank(), axis.ticks=element_blank())+
  theme(legend.position="none")+
  theme(panel.border = element_rect(fill=NA, colour = "black", linewidth=1, linetype = "solid"), panel.background = element_rect(fill="white"))

NMDSplot_occburr_elevbird # this Figure 3 in the corresponding manuscript

ggsave(filename = here("Figures", "NMDSplot_occburr_elevbird"), 
       width=7, height=5.5, units="in", dpi=300, device="png")

# Finally, we can use PERMANOVA to determine if each occupied burrow has a distinct chemical profile using:
# the compounds that are abundant on storm-petrel plumage
# adonis2() in the vegan package works in this situation because it is a balanced design where each burrow has 3 replicates
# however, we will also export the data frame and run the analysis in PRIMER to be consistent with other PERMANOVA tests in this study

set.seed(724)
burr_permanova_elevbird <- adonis2(occburr_elevbird_bray ~ Location, data = occburr_elevbird, permutations=9999)
burr_permanova_elevbird
saveRDS(burr_permanova_elevbird, here("Models", "burrow_permanova_elevatedbird.rds"))

# export data to be used in PRIMER 
writexl::write_xlsx(occburr_elevbird, here("PRIMER/Imported Files", "Occburrows_elevbirdcomps.xlsx"))

# Results from analysis in PRIMER can be viewed in: PRIMER/Exported Results/Permanova Occupied Burrows Elevated Bird Compounds.rtf