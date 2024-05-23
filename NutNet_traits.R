################################################################################
##  NutNet_traits.R: Imputing traits for NutNet.
##
##  Author: Kimberly Komatsu
##  Date created: April 17, 2024
################################################################################

library(Taxonstand)
library(WorldFlora)
library(tidyverse)


#kim's laptop
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\nutnet')


###read in data
nutnet <- read.csv('full-cover_2023-11-07.csv')%>%
  rename(cover=max_cover, genus_species=Taxon)%>%
  filter(live==1, !(genus_species %in% c('GROUND', 'OTHER LITTER', 'OTHER ARISTIDA CONTORTA (DEAD)', 'OTHER SALSOLA KALI (DEAD)', 'OTHER TRIODIA BASEDOWII (DEAD)', 'OTHER ANIMAL DROPPINGS', 'OTHER ROCK', 'OTHER ANIMAL DIGGINGS', 'OTHER WOODY OVERSTORY', 'OTHER STANDING DEAD', 'OTHER ANIMAL DIGGING', 'OTHER SOIL BIOCRUST', 'OTHER WOOD', 'OTHER SHELL', 'DEER')))%>%
  mutate(trt=as.character(ifelse(year_trt<1, 'Control', as.character(trt))))

WFO.file<-read.delim("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\WFO_Backbone\\classification.txt")


nutnetSp <- TPL(unique(nutnet$genus_species))

nutnetSpClean <- nutnetSp %>% 
  filter(!is.na(New.Species), New.Species!='sp.') %>% 
  mutate(database='NutNet', species_matched=paste(New.Genus, New.Species, sep=' ')) %>% 
  select(database, species_matched) %>% 
  unique() %>% 
  mutate(species_matched=str_to_sentence(species_matched))

GExSp <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\GEx_species_family_May2023.csv') %>% 
  select(database, species_matched) %>% 
  unique()

CoRREsp <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\TRY\\corre2trykey_2021.csv') %>% 
  select(species_matched) %>% 
  mutate(database='CoRRE') %>% 
  unique()

allSpp <- rbind(nutnetSpClean, GExSp, CoRREsp) %>% 
  pivot_wider(names_from=database, values_from=database, values_fill='not') %>% 
  mutate(nutnet_only=ifelse(GEx=='not' & CoRRE=='not', 'needs traits', 'has traits'))

nutnetFamilies <- nutnetSp %>% 
  filter(!is.na(New.Species), New.Species!='sp.') %>% 
  mutate(database='NutNet', species_matched=paste(New.Genus, New.Species, sep=' ')) %>% 
  select(database, species_matched, Family) %>% 
  unique() %>% 
  mutate(species_matched=str_to_sentence(species_matched))

traitsNeeded <- nutnet %>% 
  mutate(lifespan=str_to_lower(local_lifespan), 
         g_form=ifelse(functional_group %in% c('GRAMINOID', 'GRASS'), 'graminoid',
                     ifelse(functional_group=='LEGUME', 'forb',
                     str_to_lower(functional_group)))) %>% 
  select(genus_species, ps_path, lifespan, g_form) %>% 
  unique() %>% 
  rename(Taxon=genus_species) %>% 
  left_join(nutnetSp)  %>% 
  filter(!is.na(New.Species), New.Species!='sp.') %>% 
  mutate(database='NutNet', species_matched=paste(New.Genus, New.Species, sep=' ')) %>% 
  select(species_matched, Family, g_form, ps_path, lifespan) %>% 
  unique() %>% 
  mutate(species_matched=str_to_sentence(species_matched)) %>% 
  left_join(allSpp) %>% 
  filter(nutnet_only=='needs traits') %>% 
  mutate(alt_photopath_possible=ifelse(Family %in% c('Acanthaceae', 'Aizoaceae', 'Amaranthaceae', 'Asteraceae', 'Boraginaceae', 'Cleomaceae', 'Caryophyllaceae', 'Cyperaceae', 'Euphorbiaceae', 'Gisekiaceae', 'Hydrocharitaceae', 'Molluginaceae', 'Nyctaginaceae', 'Polygonaceae', 'Portulacaceae', 'Poaceae', 'Scrophulariaceae', 'Zygophyllaceae', 'Cactaceae', 'Crassulaceae', 'Euphorbiaceae', 'Liliaceae', 'Bromeliaceae', 'Orchidaceae'), 'possible', 'no')) %>% 
  mutate(photosynthetic_pathway=ifelse(ps_path=='NULL' & alt_photopath_possible=='no', 'C3',
                                ifelse(ps_path=='NULL' & alt_photopath_possible=='possible', 'CHECK',
                                ps_path)),
         growth_form=ifelse(g_form=='null', 'CHECK', g_form)) %>% 
  select(species_matched, Family, growth_form, photosynthetic_pathway, lifespan)

# write.csv(traitsNeeded, 'NutNet_categorical trait data_2024_to fill.csv')