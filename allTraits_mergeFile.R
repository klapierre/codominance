################################################################################
##  allTraits_mergeFile.R: Merging trait files from all three databases.
##
##  Author: Kimberly Komatsu
##  Date created: April 17, 2024
################################################################################

library(readxl)
library(PerformanceAnalytics)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data')

##### CoRRE and GEx Traits from EDI #####
correGExTraitsContinuous <- read.csv('https://portal.edirepository.org/nis/dataviewer?packageid=edi.1533.3&entityid=5ebbc389897a6a65dd0865094a8d0ffd') %>% 
  select(-source, -error_risk_overall)

correGExTraitsCategorical <- read.csv('https://portal.edirepository.org/nis/dataviewer?packageid=edi.1533.3&entityid=5ebbc389897a6a65dd0865094a8d0ffd') %>% 
  select(-source, -error_risk_overall)


##### NutNet Traits - imputed/gathered for this project #####
nutnetTraitsContinuous <- read.csv('nutnet/NutNet_continuousTraitData_imputed_20240711.csv') %>% 
  rename(family=Family) %>% 
  select(family, species, trait, trait_value)

nutnetTraitsCategorical <- read_xlsx('NutNet/NutNet_categorical_traits_2024.xlsx') %>% 
  select(family, species_matched, leaf_type, leaf_compoundness, growth_form, photosynthetic_pathway,
         lifespan, stem_support, clonal) %>% 
  rename(species=species_matched) %>% 
  pivot_longer(leaf_type:clonal, names_to='trait', values_to='trait_value')


##### rbind, pivot wider, and merge #####
continuousTraits <- rbind(correGExTraitsContinuous, nutnetTraitsContinuous) %>% 
  pivot_wider(names_from=trait, values_from=trait_value)


### start here, this is a problem becasue something is repeated and so can't pivot #####
categoricalTraits <- rbind(correGExTraitsCategorical, nutnetTraitsCategorical) %>% 
  pivot_wider(names_from=trait, values_from=trait_value) %>% 
  select(-mycorrhizal_type, -n_fixation_type)

allTraits <- continuousTraits %>% 
  full_join(categoricalTraits)


test <- categoricalTraits %>% 
  select(family, species) %>% 
  unique()
