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


##### Figuring out what proportion of cover and species we have all traits for #####
correNames <- read.csv('CoRRE\\corre2trykey_2021.csv') %>% 
  select(genus_species, species_matched) %>% 
  unique()

correAbund <- read.csv('CoRRE\\CoRRE_RawAbundance_2021.csv') %>% 
  left_join(correNames) %>% 
  rename(species2=species_matched) %>% 
  mutate(species=ifelse(is.na(species2), genus_species, species2)) %>% 
  select(site_code, project_name, community_type, block, plot_id, treatment, calendar_year, species, abundance) 

GExAbund <- read.csv('GEx\\GEx_cleaned_11June2020.csv') %>% 
  rename(site_code=site, plot_id=plot, treatment=trt, abundance=relcov, calendar_year=year) %>% 
  mutate(project_name=NA, community_type=NA) %>% 
  mutate(species=ifelse(is.na(genus_species_clean), genus_species, genus_species_clean)) %>% 
  select(site_code, project_name, community_type, block, plot_id, treatment, calendar_year, species, abundance)

###START HERE: a few experiments have multiple data points within each plot (comped multiple subplots). Need to get one value for each plot per species in some way (pick one or average)
NutNetAbund <- read.csv('NutNet\\full-cover_2023-11-07.csv') %>% 
  left_join(read_csv('NutNet\\NutNet_clean_spp_names_20240710 - Copy.csv')) %>% 
  filter(!(functional_group %in% c('BRYOPHYTE', 'NON-LIVE', 'LIVERWORT', 'LICHEN')),
         live==1) %>% 
  mutate(species2=paste(New.Genus, New.Species, sep=' ')) %>% 
  mutate(species=ifelse(species2=='NA NA', Taxon, species2)) %>% 
  rename(plot_id=plot, treatment=trt, abundance=max_cover, treatment=trt, calendar_year=year) %>% 
  mutate(project_name=NA, community_type=NA) %>% 
  select(site_code, project_name, community_type, block, plot_id, treatment, calendar_year, species, abundance)

coverData <- rbind(correAbund, GExAbund, NutNetAbund)

totCover <- coverData %>% 
  group_by(site_code, project_name, community_type, block, plot_id, calendar_year) %>% 
  summarize(tot_cover=sum(abundance)) %>% 
  ungroup()

richness <- coverData %>% 
  select(site_code, project_name, community_type, species) %>% 
  unique() %>% 
  group_by(site_code, project_name, community_type) %>% 
  summarize(richness_all=length(species)) %>% 
  ungroup()  

# relCoverAll <- coverData %>% 
#   left_join(totCover) %>% 
#   mutate(rel_cover=100*(abundance/tot_cover)) %>% 
#   ungroup() %>%
#   left_join(richness) %>% 
#   left_join(allTraits)

# relCoverWide <- relCoverAll %>% 
#   select(site_code, project_name, community_type, block, plot_id, calendar_year, 
#          treatment, species, rel_cover) %>% 
#   # group_by(site_code, project_name, community_type, block, plot_id, calendar_year,
#   #          treatment, species) %>%
#   # summarise(length=length(rel_cover)) %>%
#   # ungroup() %>%
#   pivot_wider(names_from=species, values_from=rel_cover, values_fill=0)




###proportion of species remaining if we remove NAs

codomNum <- df_grouped %>% 
  ungroup() %>% 
  select(site_code, project_name, community_type, plot_id, calendar_year, num_codominants) %>% 
  unique()

proportionSppRemain <- relCoverAll %>%
  left_join(df_grouped) %>% 
  filter(num_codominants>1, num_codominants<4) %>% 
  filter(!is.na(LDMC)) %>% 
  select(site_code, project_name, community_type, species) %>% 
  unique() %>% 
  group_by(site_code, project_name, community_type) %>% 
  summarise(richness_remaining=length(species)) %>% 
  ungroup() %>% 
  left_join(richness) %>% 
  mutate(prop_richness_remaining=richness_remaining/richness_all) %>% 
  filter(prop_richness_remaining>0.5)

hist(proportionSppRemain$prop_richness_remaining)

test <- proportionSppRemain %>% 
  select(site_code, project_name, community_type) %>% 
  unique()
#if we have a cutoff of 80% of species have traits, then we would drop 464 experiments out of 550
#cutoff of 70%, then we drop 365 experiments out of 550


### percent cover that has traits- this might not be as relevant to our dendrograms

coverWithTraits <- relCoverAll %>% 
  filter(!is.na(LDMC)) %>% 
  group_by(site_code, project_name, community_type, block, plot_id, calendar_year, richness_all) %>% 
  summarise(cover_remaining=sum(rel_cover)) %>% 
  ungroup() %>% 
  filter(cover_remaining<80)

test <- coverWithTraits %>% 
  select(site_code, project_name, community_type) %>% 
  unique()

#we have 550 experiments (site_proj_comm) and dropping any experiments that have at least one plot*year with <80% cover containing traits results in dropping 400 experiments

hist(coverWithTraits$cover_remaining)
