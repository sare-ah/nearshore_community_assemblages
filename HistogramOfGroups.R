# Summary figure of how many species for a particular species group are observed at each depth category

# Start fresh
rm(list=ls())

# Load packages
#--------------
library(tidyverse)

# Input data
#-----------

# Read in species observations
spp <- read.csv("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/SpeciesBy_matrices/AllSpeciesByDepthCategory.csv", stringsAsFactors = F)

# Read in species look-up table
lut <- read.csv("C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Data/LookupTbls/SpeciesLookUpTbl.csv", stringsAsFactors = F)

# Species composition by depth category
#--------------------------------------

# Pull out depth category & remove transdepth
spp$Depth <- str_sub(spp$TransDepth, -1)
spp$TransDepth <- NULL
summary(spp)

# Melt observations into dataframe of unique depth cat | spp
df <- spp %>% 
  pivot_longer(-Depth, names_to = "Species_Code", values_to = "count") %>% 
  dplyr::filter(count==1) %>% 
  distinct()

# Match taxonomic group to each unique observation
df <- dplyr::left_join(df, lut, by="Species_Code")

# Set Group as a factor and recode all algae groups, drop NA's (they are errors)
df$Group <- as.factor(df$Group)
df <- df %>% 
  mutate(Group = recode(Group,
                        'Brown'= 'Algae',
                        'Green'= 'Algae',
                        'Red'= 'Algae',
                        'Other'= 'Algae',
                        'Vascular plant'= 'Algae')) %>% 
  drop_na()

# Select columns of interest, calculate total counts
df <- df %>% 
  dplyr::select(Depth, count, Group) %>% 
  group_by(Group, Depth) %>% 
  mutate(Cnt = sum(count)) %>% 
  dplyr::select(Depth, Cnt, Group) %>% 
  distinct()
df

# Plot histogram - number of species per taxonomic group
p <- ggplot(df, aes(x=Depth, y=Cnt, fill=Group)) +
  geom_bar(stat="identity") +
  scale_fill_discrete()
p

# Species composition by substrate
#---------------------------------
