################################################################################
##  allTraits_mergeFile.R: Merging trait files from all three databases.
##
##  Author: Kimberly Komatsu
##  Date created: October 30, 2024
################################################################################

library(readxl)
library(PerformanceAnalytics)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data')

##### CoRRE and GEx Traits from EDI #####
correGExTraitsContinuous <- read.csv('https://portal.edirepository.org/nis/dataviewer?packageid=edi.1533.3&entityid=169fc12d10ac20b0e504f8d5ca0b8ee8') %>% 
  select(-family, -source, -imputation_error, -error_risk_overall, -error_risk_family, -error_risk_genus)

correGExTraitsCategorical <- read.csv('https://portal.edirepository.org/nis/dataviewer?packageid=edi.1533.3&entityid=5ebbc389897a6a65dd0865094a8d0ffd') %>% 
  select(-family, -source, -error_risk_overall)


##### NutNet Traits - imputed/gathered for this project #####
nutnetTraitsContinuous <- read.csv('nutnet/NutNet_continuousTraitData_imputed_20240711.csv') %>% 
  select(species, trait, trait_value) %>% 
  filter(!species %in% correGExTraitsContinuous$species)

## NutNet N-fixers ##
nutnetNfixer <- read.csv('nutnet/NutNet_species_list_N-fixers.csv') %>% 
  select(species_matched, n_fixer) %>% 
  rename(species=species_matched,
         n_fixation_type=n_fixer)

nutnetTraitsCategorical <- read_xlsx('NutNet/NutNet_categorical_traits_2024.xlsx') %>% 
  select(species_matched, leaf_type, leaf_compoundness, growth_form, photosynthetic_pathway,
         lifespan, stem_support, clonal) %>% 
  rename(species=species_matched) %>% 
  full_join(nutnetNfixer) %>% 
  pivot_longer(leaf_type:n_fixation_type, names_to='trait', values_to='trait_value') %>% 
  unique() %>% 
  filter(!species %in% correGExTraitsCategorical$species)


##### rbind, pivot wider, and merge #####
continuousTraits <- rbind(correGExTraitsContinuous, nutnetTraitsContinuous) %>% 
  pivot_wider(names_from=trait, values_from=trait_value)

categoricalTraits <- rbind(correGExTraitsCategorical, nutnetTraitsCategorical) %>% 
  pivot_wider(names_from=trait, values_from=trait_value) %>% 
  select(-mycorrhizal_type)

allTraits <- continuousTraits %>% 
  full_join(categoricalTraits)

# write.csv(allTraits, 'allTraits_CoRREGExNutNet_20241030.csv', row.names=F)