################################################################################
##  NutNet_traits.R: Imputing traits for NutNet.
##
##  Author: Kimberly Komatsu
##  Date created: April 17, 2024
################################################################################

library(Taxonstand)
library(WorldFlora)
library(data.table)
library(readxl)
library(PerformanceAnalytics)
library(scales)
library(ggpubr)
library(tidyverse)


#kim's laptop
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\nutnet')


###read in data
nutnet <- read.csv('full-cover_2023-11-07.csv')%>%
  rename(cover=max_cover, genus_species=Taxon)%>%
  filter(live==1, !(genus_species %in% c('GROUND', 'OTHER LITTER', 'OTHER ARISTIDA CONTORTA (DEAD)', 'OTHER SALSOLA KALI (DEAD)', 'OTHER TRIODIA BASEDOWII (DEAD)', 'OTHER ANIMAL DROPPINGS', 'OTHER ROCK', 'OTHER ANIMAL DIGGINGS', 'OTHER WOODY OVERSTORY', 'OTHER STANDING DEAD', 'OTHER ANIMAL DIGGING', 'OTHER SOIL BIOCRUST', 'OTHER WOOD', 'OTHER SHELL', 'DEER')))%>%
  mutate(trt=as.character(ifelse(year_trt<1, 'Control', as.character(trt))))

### only re-run if new spp needed, otherwise import from below (takes a long time to run)
# WFO.file<-read.delim("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\WFO_Backbone\\classification.txt")
# nutnetSp <- TPL(unique(nutnet$genus_species))
# write.csv(nutnetSp, 'NutNet_clean_spp_names_20240710.csv', row.names=F)

nutnetSp <- read.csv('NutNet_clean_spp_names_20240710.csv', fileEncoding="latin1")

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


##### Continuous Traits #####

#### Gather traits from databases (TRY, AusTraits, BIEN, CPTD2, TiP Leaf) ####

##### TRY data #####
dat <- fread("C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\TRY\\TRY_Traits_Downloaded_April2023.txt",sep = "\t",data.table = FALSE,stringsAsFactors = FALSE,strip.white = TRUE)

# merge NutNet with TRY
trysp <- dat %>% 
  select(AccSpeciesID, AccSpeciesName) %>% 
  rename(species_matched=AccSpeciesName) %>% 
  unique()

nutnet_key <- traitsNeeded %>% 
  right_join(trysp) %>% 
  na.omit() %>% 
  select(species_matched, AccSpeciesID, Family) %>% 
  unique()

dat2 <- dat %>%
  right_join(nutnet_key)

# selecting desired continuous traits
dat3 <- dat2 %>%
  filter(TraitID %in% c(3106,  #vegetative height
                        3109, 3110, 3114, #leaf area
                        55, #leaf dry mass
                        47, #LDMC
                        3115, 3116, 3117, #SLA
                        14, #leaf N
                        1080, 614, #SRL
                        26)) %>% #seed dry mass
  # give names to numbers for the core traits
  mutate(CleanTraitName=ifelse(TraitID==14, 'leaf_N',
                        ifelse(TraitID==26, 'seed_dry_mass', 
                        ifelse(TraitID==47, 'LDMC', 
                        ifelse(TraitID==55, 'leaf_dry_mass', 
                        ifelse(TraitID==1080, 'SRL',
                        ifelse(TraitID==3106, 'plant_height_vegetative', 
                        ifelse(TraitID==3116, 'SLA', 
                        ifelse(TraitID==3110, 'leaf_area',
                        TraitID))))))))) %>%
  filter(!is.na(StdValue)) %>% # drop observations without a trait value
  filter(is.na(OrigObsDataID)) %>% # drop known repeats in TRY 
  filter(UncertaintyName!="Range" & UncertaintyName!="Class range") # drop 1190 observations that were range estimates for leaf area, plant vegetative height and seed mass, which led to repeats in the data

# removing dead plants
health <- dat %>%
  select(DatasetID, DataID, ObsDataID, AccSpeciesID, AccSpeciesName, 
         TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, 
         StdValue, UnitName, ErrorRisk) %>%
  filter(DataID==1961) %>%
  mutate(drop=ifelse(OrigValueStr=="Dead", 1, 0)) %>%
  select(ObsDataID, drop) %>%
  unique() #list of dead plants

healthy <- dat3 %>% #merge to drop observations on dead plants
  left_join(health) %>% 
  mutate(drop=ifelse(is.na(drop), 0, drop)) %>% 
  filter(drop!=1) %>%  #no overlap in dataset
  select(-drop)

# removing trees that are not seedlings -- based on data identified specifically as either mature or seedling
treesp <- read.csv("nutnet_species_families_trees_2024.csv") #read in which species are trees

tree <- dat %>% #get list of tree observations that were made on seedlings
  select(DatasetID, DataID, ObsDataID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, 
         OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk) %>%
  filter(DataID==413) %>%
  right_join(treesp) %>%
  filter(tree.non.tree=="tree") %>%
  mutate(drop=ifelse(OrigValueStr=="seedlings"|OrigUnitStr==0|OrigValueStr=="seedling"|OrigValueStr=="Seedling (0 - 1 y)"|OrigValueStr=="seedlings, 1st year",  0, 1)) %>%
  select(ObsDataID, drop) %>%
  unique()

nontree <- dat3 %>% #merge to drop tree observations that are not seedlings
  left_join(tree) %>% 
  mutate(drop=ifelse(is.na(drop), 0, drop)) %>% 
  filter(drop!=1) %>%  #no overlap in dataset
  select(-drop)

# Removing plants that were not measured in natural conditions
setting <- dat %>% #get list of observations that were not in natural settings
  select(DatasetID, DataID, ObsDataID, AccSpeciesID, AccSpeciesName, 
         TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, 
         StdValue, UnitName, ErrorRisk) %>%
  filter(DataID==327) %>%
  mutate(drop=ifelse(OrigValueStr %in% c("Canadian High Arctic Research Station", "Control Plot", "field","Field", "Field (CG)", "Field (NE)", "field experiment", "Field Experiment", "Field plants", "forest stand", "Forest trees", "Forest understorey", "Fully open overstory 90 days, seedling", "Fully open overstory, seedling","Fully sunlit - Natural environment","High desert", "in situ", "In situ", "La Selva Biological Station", "meadows (M) and pastures (P) on south east to south west exposed slopes", "Montane meadow", "Mosses in forest", "nat env", "natural", "Natural", "natural-environment", "natural env", "natural envireonment", "natural enviroment", "Natural Enviroment", "natural environment", "Natural environment", "Natural Environment", "natural environment, high regional N and S deposition","natural environment, no warming, preccipitation ambient", "natural environment, sun exposed", "Natural Envrionment", "Natural Forest", "natural forest environment", "natural vegetation", "Natural Vegetation", "Natural Vegetation", "natural vegetation, but not top canopy", "natural wetland environment", "natural wetlands (field conditions)", "Natural/C", "natural_environment", "none", "None", "North facing slope", "Shade - Natural environment","South facing slope", "Trees in field"), 0, 1)) %>%
  select(ObsDataID, drop) %>%
  unique()

natural <- dat3 %>% # merge to drop observations in non-natural settings
  left_join(setting) %>% 
  mutate(drop=ifelse(is.na(drop), 0, drop)) %>% 
  filter(drop!=1) %>%  #no overlap in dataset
  select(-drop)

# Drop traits for trees -- based on whether or not the species is a tree (for all that were not designated as seedling, above)
splist <- nutnet_key %>% 
  select(species_matched) %>%
  unique() %>%
  left_join(treesp) %>%
  unique()

cont_traits <- dat3 %>%
  left_join(splist) %>%
  mutate(remove=ifelse(tree.non.tree=="tree", 1, 0)) %>%
  filter(remove==0) %>%
  select(-tree.non.tree, -remove) # drops 6694 observations

# add taxonomic information for each species
cont_traits2 <- cont_traits %>%
  select(DatasetID, ObservationID, Family, species_matched, CleanTraitName, StdValue, ErrorRisk, Reference, UnitName, OriglName, TraitID) %>%
  separate(remove = F, species_matched, into = c("genus", "species"), sep=" ") %>%
  select(-species)

# removing trait outliers based on TRY's Error Risk designation
cont_traits3a <- cont_traits2 %>%
  select(DatasetID, ObservationID, Family, genus, species_matched, CleanTraitName, StdValue, ErrorRisk, Reference, UnitName, OriglName, TraitID) %>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk)) %>%
  filter(ErrorRisk2<3) %>% #removes all observations that are greater than 3 sd from full database mean: drops 3479 observations
  select(-ErrorRisk, -ErrorRisk2) %>% 
  filter(StdValue>0) #removing negative and 0 values (drops 8 observations)

d415 <- cont_traits3a %>% 
  filter(DatasetID==415) %>% 
  group_by(DatasetID, species_matched, CleanTraitName, Family, genus, Reference)%>%
  summarise(StdValue=mean(StdValue)) %>% 
  ungroup() %>% 
  mutate(ObservationID=row_number())

TRYtraits <- cont_traits3a %>%
  filter(!(DatasetID %in% c(415))) %>%
  bind_rows(d415) %>% 
  mutate(DatabaseID="TRY") %>% 
  select(DatabaseID, DatasetID, ObservationID, Family, species_matched, genus, CleanTraitName, StdValue, Reference) %>% 
  filter(StdValue>0)

# finds repeats; everything left is probably real, just measurement imprecision
# repeats <- cont_traits4 %>% 
#   group_by(species_matched, CleanTraitName, StdValue) %>%
#   summarize(n=length(StdValue)) %>%
#   ungroup() %>% 
#   filter(n>2)

# write.csv(TRYtraits, "NutNet_TRY traits_20240711.csv", row.names=F)


##### BIEN traits #####

#gather data from library
nutnetSpp <- nutnet_key %>% 
  select(-AccSpeciesID)

library(BIEN)

sp.vector <- unique(nutnet_key$species_matched)

bienData <- BIEN_trait_species(species=sp.vector) %>% 
  rename(species_matched=scrubbed_species_binomial) %>%  
  right_join(nutnetSpp) %>% 
  # subset to data that we want
  filter(trait_name %in% c('leaf area', 'leaf area per dry mass', 'leaf dry mass', 
                           'leaf dry mass per leaf fresh mass', 'seed mass',
                           'leaf nitrogen content per leaf dry mass')) %>% 
  mutate(trait_value=as.numeric(trait_value)) %>% 
  # standardize units to fit TRY
  mutate(clean_trait_value=ifelse(trait_name=='leaf dry mass per leaf fresh mass', trait_value/1000, #LDMC (BIEN mg/g   TRY g/g)
                           ifelse(trait_name=='leaf area per leaf dry mass', trait_value*1000, #SLA (BIEN m2/kg   TRY mm2/g)
                           ifelse(trait_name=='leaf dry mass', trait_value*1000, #leaf dry mass (BIEN g   TRY mg)
                           trait_value)))) %>% 
  # remove data that was not from a naturally growing plant
  filter(method!='laboratory/greenhouse/garden experiment',
         trait_value!=0) %>% 
  # Problem: A few datasets have lots of repeated data for some traits*species.
  # Solution: For each species, find if there is repeated data for all traits collected on an individual. 
  # Where this occurs, keep the lowest ObservationID.
  select(species_matched, trait_name, project_pi, id, clean_trait_value) %>% 
  pivot_wider(names_from=trait_name, values_from=clean_trait_value, names_prefix = "d__") %>% 
  group_by_at(vars(!id)) %>% 
  mutate(n=length(species_matched), obid2=min(id)) %>% 
  ungroup() %>% 
  select(-id) %>% 
  unique() %>% 
  pivot_longer(3:7, names_to = "CleanTraitName1", values_to = "StdValue") %>% 
  separate(CleanTraitName1, into = c("prefix", "CleanTraitName"), "__") %>% 
  select( -prefix, -n) %>% 
  na.omit() %>% 
  rename(ObservationID=obid2)

# checking for duplicate data
test <- bienData %>% 
  group_by(species_matched, CleanTraitName, StdValue) %>% 
  summarize(n=length(StdValue)) %>% 
  ungroup() %>% 
  filter(n>2)

# change BIEN trait names to fit TRY trait names
bienData$CleanTraitName <- recode(bienData$CleanTraitName, 
                                    'leaf area'='leaf_area',
                                    'leaf area per dry mass'='SLA',
                                    'leaf dry mass'='leaf_dry_mass',
                                    'leaf dry mass per leaf fresh mass'='LDMC',
                                    'leaf nitrogen content per leaf dry mass'='leaf_N',
                                    'seed mass'='seed_dry_mass')


# unify with other dataset columns
BIENtraits <- bienData %>% 
  mutate(DatabaseID='BIEN') %>% 
  rename(DatasetID=project_pi) %>% 
  select(DatabaseID, DatasetID, ObservationID, species_matched, CleanTraitName, StdValue) %>% 
  # Make a genus column
  separate(species_matched, into = c("genus","species"), sep=" ", remove=FALSE) %>% 
  select(-species) %>% 
  # Filter outliers
  mutate(drop=ifelse(DatasetID %in% c('Abakumova M', 'Liu Y', 'Osborne CP') & CleanTraitName=='leaf_area', 1, #these studies did something other than leaf area (e.g., total leaf area for the whole plant)
                     ifelse(DatasetID %in% c('Schmid B') & CleanTraitName=='leaf_dry_mass', 1, 0))) %>% #this study did something other than leaf dry mass (e.g., plant mass)
  filter(drop==0) %>% 
  select(-drop) %>% 
  left_join(nutnetSpp) %>%
  mutate(Reference=DatasetID) %>% 
  select(DatabaseID, DatasetID, ObservationID, Family, species_matched, genus, CleanTraitName, StdValue, Reference) %>% 
  filter(StdValue>0)

# write.csv(BIENtraits, 'NutNet_BIEN traits_20240711.csv', row.names=F)


##### AusTraits #####
library(austraits)

austraits <- load_austraits(version = "6.0.0", path = "austraits")

traits <- summarise_austraits(austraits, "trait_name")

#doesn't have seed number, stem specific density, rooting depth
data <- extract_trait(austraits, c('leaf_area',
                                   'leaf_dry_mass', 
                                   'leaf_dry_matter_content', 
                                   'leaf_mass_per_area', #need to inverse this
                                   'leaf_N_per_dry_mass', 
                                   'plant_height', 
                                   'root_specific_root_length', 
                                   'seed_dry_mass'))


traitData <- data$traits %>%
  mutate(DatabaseID='AusTraits') %>%
  rename(DatasetID=dataset_id,
         ObservationID=observation_id,
         species_matched=taxon_name,
         StdValue=value) %>% 
  mutate(StdValue=ifelse(trait_name=='leaf_mass_per_area', (1/StdValue)*1000, 
                  ifelse(trait_name=='root_specific_root_length', (StdValue*100), StdValue)),
         trait_name=ifelse(trait_name=='leaf_mass_per_area', 'specific_leaf_area', trait_name))

species <- nutnetSpp %>%  #species names are standardized
  left_join(treesp) %>% 
  filter(tree.non.tree != "tree") %>% #Remove trees
  separate(species_matched, into=c('genus', 'species', 'subspp'), sep=' ') %>% 
  filter(species!='sp.') %>% 
  unite(col='species_matched', genus:species, sep=' ', remove=T) %>% 
  select(species_matched) %>% 
  unique()

AusTraits <- traitData %>%
  right_join(species) %>%
  group_by(DatabaseID, DatasetID, ObservationID, species_matched, trait_name) %>% 
  summarize(StdValue=max(StdValue)) %>% 
  ungroup() %>% 
  mutate(drop=ifelse(trait_name=='plant_height' & StdValue>40, 1,
                     ifelse(trait_name=='seed_dry_mass' & StdValue>600, 1, 0))) %>% 
  filter(drop==0) %>% 
  select(-drop)

AusTraits$CleanTraitName <- recode(AusTraits$trait_name, 
                                        'leaf_dry_matter_content'='LDMC', 
                                        'specific_leaf_area'='SLA', 
                                        'leaf_N_per_dry_mass'='leaf_N', 
                                        'plant_height'='plant_height_vegetative', 
                                        'root_specific_root_length'='SRL') 
AusTraits <- AusTraits %>%
  left_join(nutnetSpp) %>%
  mutate(species_matched2=species_matched) %>%
  separate(species_matched2, into=c('genus', 'species')) %>%
  mutate(DatabaseID='AusTraits',
         Reference=DatasetID) %>%
  select(DatabaseID, DatasetID, ObservationID, Family, species_matched, genus, CleanTraitName, StdValue, Reference) %>% 
  filter(StdValue>0)

# write.csv(AusTraits, 'NutNet_AusTraits_20240711.csv', row.names=F)

##### TiP Leaf #####
tipTraits <- read_xlsx('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\TiP_leaf\\The TiP-Leaf dataset.xlsx', sheet='plant traits') %>% 
  dplyr::rename(species_matched=Species) %>%
  filter(species_matched!='/') %>% 
  mutate(ObservationID=row_number(.)) %>% 
  mutate(DatasetID='1', DatabaseID='TIPleaf') %>% 
  left_join(nutnetSpp) %>% 
  filter(!is.na(Family)) %>% 
  mutate(LCC=as.numeric(ifelse(LCC=='/', NA, LCC)),
         LNC=as.numeric(ifelse(LNC=='/', NA, LNC)),
         LPC=as.numeric(ifelse(LPC=='/', NA, LPC)),
         SLA=SLA/10) %>% #unit conversion to TRY standards: cm2/g to mm2/mg 
  select(DatabaseID, DatasetID, ObservationID, species_matched, DW, LDMC, LA, SLA, LNC) %>% 
  pivot_longer(DW:LNC, names_to='trait_name', values_to='StdValue') %>% 
  unique()

tipTraits$CleanTraitName <- recode(tipTraits$trait_name, 
                             'DW'='leaf_dry_mass',
                             'LA'='leaf_area',
                             'LNC'='leaf_N')

tipTraits <- tipTraits %>% 
  separate(col=species_matched, into=c('genus', 'species'), sep=' ', remove=F) %>% 
  left_join(nutnetSpp) %>% 
  mutate(Reference='TipLeaf Database') %>% 
  select(DatabaseID, DatasetID, ObservationID, Family, genus, species_matched, CleanTraitName, StdValue, Reference) %>% 
  filter(StdValue>0)

# write.csv(tipTraits, 'NutNet_TiP Leaf traits_20240711.csv', row.names=F)

##### China Plant Trait Database 2 #####
spList <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\ChinaPlant2\\Species translations.csv') %>% 
  unite(col='species_matched', ACCEPTED.GENUS:ACCEPTED.SPECIES, sep=' ') %>% 
  select(species_matched, Site.ID, SAMPLE.ID) %>% 
  left_join(nutnetSpp)

chem <- read.csv("C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\ChinaPlant2\\Chemical traits.csv") %>% 
  filter(flagged=="") %>% 
  mutate(leaf_area=Average.LA*1000000) %>% #unit conversion to TRY standards: m2 to mm2
  mutate(LDMC=LDMC/1000) %>% #unit conversion to TRY standards: mg/g to g/g
  select(-LMA,-Narea, -Parea, -Karea, -d13C.12C, -d15N.14N, -flagged, -Average.LA) %>% 
  rename(leaf_C=Cmass, 
         leaf_N=Nmass,
         leaf_P=Pmass, 
         leaf_K=Kmass) %>% 
  pivot_longer(SLA:leaf_area, names_to="CleanTraitName", values_to="StdValue") %>% 
  right_join(spList) %>% 
  na.omit()

photo <- read.csv("C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\ChinaPlant2\\Photosynthetic traits.csv") %>% 
  filter(flagged=="") %>% 
  select(SAMPLE.ID, Vcmax, Jmax) %>% 
  rename(Vc_max=Vcmax,
         J_max=Jmax) %>% 
  pivot_longer(Vc_max:J_max, names_to='CleanTraitName', values_to='StdValue') %>% 
  right_join(spList) %>% 
  na.omit()

#bind together
CPTDtraits <- rbind(chem, photo) %>% 
  mutate(DatabaseID='CPTD2') %>% 
  rename(DatasetID=Site.ID,
         ObservationID=SAMPLE.ID) %>% 
  filter(CleanTraitName %in% c('LDMC', 'leaf_area', 'leaf_N', 'SLA')) %>% 
  separate(col=species_matched, into=c('genus', 'species'), sep=' ', remove=F) %>% 
  mutate(Reference='China Plant Trait Database 2') %>% 
  select(DatabaseID, DatasetID, ObservationID, Family, genus, species_matched, CleanTraitName, StdValue, Reference) %>% 
  filter(StdValue>0)


##### Bind all trait data #####
allTraits <- rbind(TRYtraits, AusTraits, BIENtraits, tipTraits, CPTDtraits) %>% 
  mutate(ReferenceID=paste(DatabaseID, DatasetID, sep='_')) %>% 
  select(DatabaseID, DatasetID, ObservationID, Family, genus, species_matched, 
         CleanTraitName, StdValue) %>% 
  group_by(DatabaseID, DatasetID, ObservationID, Family, genus, species_matched, 
           CleanTraitName) %>% 
  summarise(StdValue=mean(StdValue)) %>% 
  ungroup()

allTraits_wide <- allTraits %>% 
  pivot_wider(names_from = CleanTraitName, values_from = StdValue, values_fill=NA)

ntraits <- length(unique(allTraits$CleanTraitName))
miss <- sum(is.na(allTraits_wide))
total <- nrow(allTraits_wide)*ntraits
miss/total*100
#missing 89.03% of data

spnum <- length(unique(allTraits_wide$species_matched))
famnum <- length(unique(allTraits_wide$Family))

label <- allTraits %>%
  group_by(CleanTraitName, DatabaseID) %>%
  summarise(length=length(StdValue)) %>%
  ungroup() %>%
  group_by(CleanTraitName) %>%
  mutate(length2=sum(length)) %>%
  ungroup() %>%
  pivot_longer(cols=length:length2, names_to='name', values_to='length') %>%
  mutate(DatabaseID=ifelse(name=='length2', 'total', DatabaseID)) %>%
  unique() %>%
  mutate(percent=round((length/253224)*100, 1)) %>% 
  mutate(CleanTraitName2=ifelse(CleanTraitName==3109, 'Leaf Area (leaflet, -petiole)',
                                ifelse(CleanTraitName==3114, 'Leaf Area (undefined, undefined)',
                                ifelse(CleanTraitName=='leaf_area', 'Leaf Area (leaf, +petiole)',
                                ifelse(CleanTraitName==3115, 'Specific Leaf Area (-petiole)',
                                ifelse(CleanTraitName==3117, 'Specific Leaf Area (undefined)',
                                ifelse(CleanTraitName=='SLA', 'Specific Leaf Area (+petiole)', 
                                ifelse(CleanTraitName=='SRL', 'Specific Root Length (all root)',
                                ifelse(CleanTraitName==614, 'Specific Root Length (fine root)', 
                                ifelse(CleanTraitName=='leaf_N', 'Leaf N Content',
                                ifelse(CleanTraitName=='plant_height_vegetative', 'Plant Vegetative Height',
                                ifelse(CleanTraitName=='seed_dry_mass', 'Seed Dry Mass',
                                ifelse(CleanTraitName=='leaf_dry_mass', 'Leaf Dry Mass',
                                ifelse(CleanTraitName=='LDMC', 'Leaf Dry Matter Content',
                                CleanTraitName))))))))))))))

label$CleanTraitName2 = factor(label$CleanTraitName2, levels=c('Leaf Area (leaf, +petiole)', 'Leaf Area (leaflet, -petiole)', 'Leaf Area (undefined, undefined)', 'Leaf Dry Mass', 'Leaf Dry Matter Content', 'Specific Leaf Area (+petiole)', 'Specific Leaf Area (-petiole)', 'Specific Leaf Area (undefined)', 'Leaf N Content', 'Plant Vegetative Height', 'Specific Root Length (all root)', 'Specific Root Length (fine root)', 'Seed Dry Mass'))

# How many observations do we have for each trait across our database?
ggplot(data=label, aes(x=DatabaseID, y=length, label=round(percent,1), fill=DatabaseID)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept=253224*.2) + # 20% of observations missing any given trait
  geom_hline(yintercept=253224*.1, color='red') + # 10% of observations missing any given trait
  geom_text(vjust = -0.25, size=6) +
  facet_wrap(~CleanTraitName2, ncol=5, labeller=label_wrap_gen(width=25)) +
  scale_y_continuous() +
  scale_x_discrete(breaks=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY", "total"),
                   limits=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY", "total"),
                   labels=c("Au", "BN", "C2", "TP", "TY", 'all')) +
  scale_fill_manual(values=c('#4E3686', '#5DA4D9', '#80D87F', '#FED23F','darkgrey', '#EE724C'))+
  theme(strip.text.x = element_text(size = 18),
        axis.title.x=element_text(size=24, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=22),
        axis.title.y=element_text(size=24, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=22),
        legend.position='none') +
  ylab('Number of Observations') + xlab('Database ID')
# ggsave('x.png', width=17, height=19, units='in', dpi=300, bg='white')

# Are there any outlier datasets for each trait?
ggplot(data=allTraits, aes(x=DatabaseID, y=StdValue)) +
  scale_y_log10() + # note log axis!
  geom_jitter(aes(color=DatabaseID)) +
  geom_boxplot(color='black', alpha=0) +
  facet_wrap(~CleanTraitName, scales='free_y', ncol=4) +
  scale_x_discrete(breaks=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY"),
                   labels=c("A", "B", "C", "TIP", "TY")) +
  scale_color_manual(values=c('#4E3686', '#5DA4D9', '#80D87F', '#FED23F', '#EE724C')) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position='top') 
# ggsave('x.png', width=7.5, height=10, units='in', dpi=300, bg='white')


# Transpose to wide format for gap filling.
talltraits <- allTraits %>% 
  group_by(DatabaseID, DatasetID, ObservationID, Family, genus, species_matched) %>%
  pivot_wider(names_from=CleanTraitName, values_from=StdValue, values_fill=NA) %>% 
  ungroup()

# write.csv(allTraits, 'nutnet_trait database_combo_continuous_20240710_long.csv', row.names = F)

# write.csv(talltraits, 'nutnet_trait database_combo_continuous_20240710.csv', row.names = F)



##### Impute Traits #####

library(BHPMF)
library(plyr)
library(abind)
library(mice)

traits <- read.table('nutnet_trait database_combo_continuous_20240710.csv', row.names=NULL, sep=",", header=T) %>% 
  mutate(family=Family) %>% 
  select(-Family)

# traits <- traits[1:100000,]

#remove trait values with > 4 SD:
spp <- unique(traits$species_matched) #get vector with species names

out<-NULL
for(i in 1:length(spp)) { #loop for each species
  print(i/length(spp))
  sub <- traits[traits$species_matched %in% spp[i],]
  for(j in 7:ncol(traits)) { #loop for each trait (column)
    mean.sub <- mean(sub[,j], na.rm=T)
    sd.sub <- sd(sub[,j], na.rm=T)
    
    lim_up <- mean.sub + 4*sd.sub
    lim_dn <- mean.sub - 4*sd.sub
    
    sub[,j][sub[,j] > lim_up] <- NA
    sub[,j][sub[,j] < lim_dn] <- NA
  }
  out <- rbind(out, sub)
}

#create hierarchy file:
hierarchy.info <- subset(traits, select = c(ObservationID, species_matched, genus, family))
names(hierarchy.info) <- c("plant_id","species", "genus", "family")
hierarchy.info$plant_id <- 1:nrow(hierarchy.info)

#some genera are assigned to different families. Need to be unified:
hierarchy.info$family[hierarchy.info$genus=="Lancea"] <- "Mazaceae"
hierarchy.info$family[hierarchy.info$genus=="Toxicoscordion"] <- "Melanthiaceae"
hierarchy.info$family[hierarchy.info$genus=="Heliotropium"] <- "Boraginaceae"
hierarchy.info$family[hierarchy.info$genus=="Phacelia"] <- "Boraginaceae"
hierarchy.info$family[hierarchy.info$genus=="Pholistoma"] <- "Boraginaceae"

# test <- hierarchy.info %>%
#   select(family, genus) %>% 
#   unique(.) %>% 
#   group_by(genus) %>%
#   summarize(length(family)) %>%
#   ungroup()

#create trait info file:
trait.info <- as.data.frame(subset(traits, select = -c(family, genus, species_matched, ObservationID,
                                                       DatabaseID, DatasetID)))

#check if both datasets are equal
nrow(hierarchy.info) == nrow(trait.info)



##### z % log transform #####
back_trans_pars <- list()
rm_col <- c()
for(i in 1:ncol(trait.info)){
  x <- trait.info[,i] # goes through the columns
  min_x <- min(x,na.rm = T) # takes the min of each column
  if(min_x < 0.00000000001){
    x <- x - min_x + 1 # make this optional if min x is neg
  }
  logx <- log10(x)
  mlogx <- mean(logx, na.rm = T)
  slogx <- sd(logx, na.rm = T)
  x <- (logx - mlogx)/slogx # Z transformation
  back_trans_pars[[i]] <- list(min_x = min_x,
                               mlogx = mlogx,
                               slogx = slogx)
  trait.info[,i] <- x
}

# write.table(back_trans_pars, "imputation_20240711\\back_trans_pars.csv")


##### gap-filling #####
#set-directory
tmp.dir <- dirname("imputation_20240711\\tmp")

#set parameters
smpl <- 900:1000
fold <- c(rep(10:20, 8), 10, 11)

#set number of iterations:
repe <- 90 #should be 90

for(i in 1:repe) { #loop for each trait (column)
  set.seed(123)
  GapFilling(as.matrix(trait.info), hierarchy.info,
             num.samples = smpl[i], num.folds.tuning=fold[i], burn=187,
             mean.gap.filled.output.path = paste0(tmp.dir,"/mean_gap_filled_",i,".txt"),
             std.gap.filled.output.path = paste0(tmp.dir,"/std_gap_filled_",i,".txt"),
             tmp.dir = tmp.dir, verbose=F)
}


##### load imputed traits and clean-up table #####
mean.trait <- list()
std.trait <- list()
for(i in 1:repe) { #loop for each trait (column)
  print(i)
  trt <- read.table(paste0("imputation_20240711\\mean_gap_filled_",i,".txt"), row.names=NULL, header=T)
  std <- read.table(paste0("imputation_20240711\\std_gap_filled_",i,".txt"), row.names=NULL, header=T)
  
  #Return to NA those values with SD > 1:
  for(j in 1:ncol(trt)) {
    trt[,j][std[,j]>1] <- NA
  }
  
  #Return to NA values > 1.5*max observed trait:
  for(j in 1:ncol(trt)) {
    maxt <- max(trait.info[,j], na.rm=T)
    trt[,j][trt[,j] > (maxt*1.5)] <- NA
  }
  
  mean.trait[[i]] <- trt
  std.trait[[i]] <- std
}


#### get mean across all means ####
mean.trait <- abind(mean.trait, along=3)
mean.trait <- apply(mean.trait, c(1,2), mean, na.rm=T)
mean.trait[is.nan(mean.trait)] <- NA

#data for back transforming output
back <- read.table("imputation_20240711\\back_trans_pars.csv")

#don't replace original values:
trait.info.noreplacement <- as.data.frame(mean.trait)

o <- 1 #to select the appropriate columns:
for(i in 1:ncol(trait.info.noreplacement)){
  
  #recover values:
  min_x <- back[1,o]
  mlogx <- back[1,o+1]
  slogx <- back[1,o+2]
  
  #back transform:
  x <- trait.info.noreplacement[,i] # goes through the columns
  logx <- (x*slogx) + mlogx
  b <- 10^logx
  
  #for negative values
  if(min_x < 0.00000000001){
    b <- b + min_x - 1 # make this optional if min x is neg
  }
  
  trait.info.noreplacement[,i] <- b
  o <- o+3
}

#save output
# write.csv(trait.info.noreplacement, "imputation_20240711\\imputed_traits.csv", row.names=F)


#### get mean across all std ####
std.trait <- abind(std.trait, along=3)
std.trait <- apply(std.trait, c(1,2), mean, na.rm=T)
std.trait[is.nan(std.trait)] <- NA

#don't replace original values:
trait.std.noreplacement <- as.data.frame(std.trait)

o <- 1 #to select the appropriate columns:
for(i in 1:ncol(trait.std.noreplacement)){
  
  #recover values:
  min_x <- back[1,o]
  mlogx <- back[1,o+1]
  slogx <- back[1,o+2]
  
  #back transform:
  x <- trait.std.noreplacement[,i] # goes through the columns
  logx <- (x*slogx) + mlogx
  b <- 10^logx
  
  #for negative values
  if(min_x < 0.00000000001){
    b <- b + min_x - 1 # make this optional if min x is neg
  }
  
  trait.std.noreplacement[,i] <- b
  o <- o+3
}

#save output
# write.csv(trait.std.noreplacement, "imputation_20240711\\imputed_traits_std.csv", row.names=F)



##### Impute missing values with "mice" #####
trait.info.mice <- complete(mice(trait.info.noreplacement, method="cart"), action = "long")
trait.info.mice.mean <- aggregate(. ~ .id, data = trait.info.mice[, -1], FUN = mean) #mean values: the final output
trait.info.mice.sd <- aggregate(. ~ .id, data = trait.info.mice[, -1], FUN = sd) #SDs per observation

# write.csv(trait.info.mice.mean, "imputation_20240711\\imputed_traits_mice.csv", row.names=F)
# write.csv(trait.info.mice.sd, "imputation_20240711\\imputed_traits_mice_std.csv", row.names=F)

#clean-up:
# rm(list = ls())



##### Clean-up imputed traits #####

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=26),
             axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=26),
             plot.title = element_text(size=54, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=30))


#### Categorical trait data ####
categoricalTraits <- read.csv("XYZ.csv") %>% 
  filter(species_matched!='') %>% 
  dplyr::select(family, species_matched, leaf_type, leaf_compoundness, stem_support, growth_form, photosynthetic_pathway, 
                lifespan,  clonal, mycorrhizal_type, n_fixation, rhizobial, actinorhizal, leaf_type_source, 
                leaf_compoundness_source, growth_form_source, photosynthetic_pathway_source, lifespan_source, 
                stem_support_source, clonal_source, mycorrhizal_source, n_fixation_source) %>%
  mutate(photosynthetic_pathway = replace(photosynthetic_pathway, grep("possible", photosynthetic_pathway), "uncertain")) %>%
  # mutate(clonal = replace(clonal, clonal=="uncertain", NA)) %>%
  mutate(mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type %in% c("arbuscular", "facultative_AM"), "AM"),
         mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type %in% c("double_AM_EcM", "EcM-AM", "facultative_AM_EcM", 
                                                                              "NM-AM", "NM-AM, rarely EcM",
                                                                              "species-specific: AM or rarely EcM-AM or AM"), "multiple"),
         mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="ecto", "EcM"),
         mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="ericaceous", "ErM"),
         mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="orchidaceous", "OM"),
         mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="Thysanothus", "TM"),
         mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="NM", "none"),
         mycorrhizal_type = replace(mycorrhizal_type, is.na(mycorrhizal_type), "uncertain")) %>%
  # mutate(lifespan = replace(lifespan, lifespan=="uncertain", NA)) %>%
  mutate(n_fixation_type=ifelse(rhizobial=='yes', 'rhizobial',
                         ifelse(actinorhizal=='yes', 'actinorhizal', 'none'))) %>% 
  filter(lifespan != "moss") %>% 
  select(-n_fixation, -rhizobial, -actinorhizal)

categorical_TRY <- categoricalTraits %>% 
  select(leaf_type_source, leaf_compoundness_source, growth_form_source, photosynthetic_pathway_source, lifespan_source, stem_support_source, clonal_source, mycorrhizal_source, n_fixation_source) %>% 
  pivot_longer(leaf_type_source:n_fixation_source, names_to='trait', values_to='source') %>% 
  filter(grepl("TRY", source))

categoricalTraitsFamilies <- rbind(categoricalTraits, categoricalTraitsGEx) %>% 
  filter(lifespan != "moss") %>% 
  filter(species_matched!='') %>% 
  select(family) %>% 
  unique()

categoricalTraitsError <-  rbind(categoricalTraits, categoricalTraitsGEx) %>% 
  filter(lifespan != "moss") %>% 
  filter(species_matched!='') %>% 
  select(species_matched, growth_form_error, photosynthetic_pathway_error, lifespan_error, stem_support_error, clonal_error) %>% 
  filter(growth_form_error!='')


#### Pie Charts for each categorical trait ####
# leaf type
leafType <- categoricalTraits %>% 
  group_by(leaf_type) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

leafTypeFig <- ggplot(leafType, aes(x="", y=proportion, fill=leaf_type)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF'))  +
  ggtitle('(d) Leaf Type') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# leaf compoundness
leafCompoundness <- categoricalTraits %>% 
  group_by(leaf_compoundness) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

leafCompoundnessFig <- ggplot(leafCompoundness, aes(x="", y=proportion, fill=leaf_compoundness)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF'))  +
  ggtitle('(e) Leaf Compoundness') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# stem support
stemSupport <- categoricalTraits %>% 
  group_by(stem_support) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

stemSupportFig <- ggplot(stemSupport, aes(x="", y=proportion, fill=stem_support)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF'))  +
  ggtitle('(f) Stem Support') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# growth form
growthForm <- categoricalTraits %>% 
  group_by(growth_form) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

growthFormFig <- ggplot(growthForm, aes(x="", y=proportion, fill=growth_form)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF')) +
  ggtitle('(a) Growth Form') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# photosynthetic pathway
photosyntheticPathway <- categoricalTraits %>% 
  group_by(photosynthetic_pathway) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

photoPathFig <- ggplot(photosyntheticPathway, aes(x="", y=proportion, fill=photosynthetic_pathway)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF'))  +
  ggtitle('(g) Photosynthetic Path') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# lifespan
lifespan <- categoricalTraits %>% 
  group_by(lifespan) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

lifespanFig <- ggplot(lifespan, aes(x="", y=proportion, fill=lifespan)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF'))  +
  ggtitle('(b) Lifespan') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# clonal
clonal <- categoricalTraits %>% 
  group_by(clonal) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

clonalFig <- ggplot(clonal, aes(x="", y=proportion, fill=clonal)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF')) +
  ggtitle('(c) Clonality') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# mycorrhizal type
mycorrhizalType <- categoricalTraits %>% 
  group_by(mycorrhizal_type) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

mycorrFig <- ggplot(mycorrhizalType, aes(x="", y=proportion, fill=mycorrhizal_type)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF')) +
  ggtitle('(h) Mycorrhizal Type') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

# n fixation type
nFixationType <- categoricalTraits %>% 
  group_by(n_fixation_type) %>%
  count() %>% 
  ungroup() %>% 
  mutate(proportion = round((n/sum(n)), digits=3)) %>% 
  arrange(proportion) %>%
  mutate(labels=scales::percent(proportion))

nFixFig <- ggplot(nFixationType, aes(x="", y=proportion, fill=n_fixation_type)) +
  geom_col() +
  coord_polar(theta="y") +
  scale_fill_manual(values=c('#7DCBBB', '#FFFFA4', '#B0AAD1', '#F7695F', '#6EA1C9', '#FBA550', '#A5DA56', '#AD68AF'))  +
  ggtitle('(i) N Fixation Type') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        plot.title = element_text(vjust = 0.5),
        legend.position = 'none')

# ggsave('xyz.png', width=8, height=8, units='in', dpi=300, bg='white')

#grouped figure
ggarrange(growthFormFig, lifespanFig, clonalFig,
          leafTypeFig, leafCompoundnessFig, stemSupportFig,
          photoPathFig, mycorrFig, nFixFig,
          ncol = 3, nrow = 3)

# ggsave('xyz2.png', width=26, height=26, units='in', dpi=300, bg='white')


#### Continuous traits ####

# Import imputed trait data and bind on species information
## this is trait data without replacement (all imputed)
imputedMean <- read.csv("imputation_20240711\\imputed_traits_mice.csv") %>%
  bind_cols(read.csv('nutnet_trait database_combo_continuous_20240710.csv')[,c('DatabaseID', 'DatasetID', 'ObservationID', 'Family', 'genus', 'species_matched')])  %>% 
  pivot_longer(names_to='trait', values_to='imputed_value', plant_height_vegetative:X3114) %>% 
  select(-.id)

# Import std deviations from BHPMF imputation
imputedStdBHPMF <- read.csv("imputation_20240711\\imputed_traits_std.csv") %>%
  bind_cols(read.csv('nutnet_trait database_combo_continuous_20240710.csv')[,c('DatabaseID', 'DatasetID', 'ObservationID', 'Family', 'genus', 'species_matched')]) %>% 
  pivot_longer(names_to='trait', values_to='std_BHPMF', plant_height_vegetative:X3114)

# Import std deviations from mice imputation
imputedStdMice <- read.csv("imputation_20240711\\imputed_traits_mice_std.csv") %>%
  bind_cols(read.csv('nutnet_trait database_combo_continuous_20240710.csv')[,c('DatabaseID', 'DatasetID', 'ObservationID', 'Family', 'genus', 'species_matched')])  %>% 
  pivot_longer(names_to='trait', values_to='std_mice', plant_height_vegetative:X3114) %>% 
  select(-.id)

# Replace std from BHPMF where value was set to NA because outside of error bounds with mice std
imputedStd <- left_join(imputedStdBHPMF, imputedStdMice) %>% 
  mutate(std=ifelse(std_mice>0, std_mice, std_BHPMF),
         imputation_method=ifelse(std_mice>0, 'MICE', 'BHPMF'))

# Merge means and std together
imputedRaw <- left_join(imputedMean, imputedStd) %>%   
  dplyr::select(-std_mice, -std_BHPMF) %>% 
  filter(trait %in% c('LDMC', 'leaf_area', 'leaf_dry_mass', 'leaf_N', 'plant_height_vegetative', 'seed_dry_mass', 'SLA', 'SRL'))


# Calculate averages for each species
meanContinuous <- imputedRaw %>% 
  group_by(Family, species_matched, trait) %>% 
  summarize_at(.vars=c('imputed_value', 'std'),
               .funs=list(mean=mean, sd=sd),
               na.rm=T) %>% 
  ungroup()

speciesCount <- meanContinuous %>% 
  select(Family, species_matched) %>% 
  unique() %>% 
  group_by(Family) %>% 
  summarize(num_species=length(Family)) %>% 
  ungroup() #117 families, 1147 species


#### Clean imputed continuous trait data ####
# Checked to ensure no negative values (confirmed that there are none)

transformed <- imputedRaw %>% 
  group_by(trait) %>% 
  mutate(log=log10(imputed_value)) %>% 
  ungroup() 

# with(subset(transformed, trait=='LDMC'), hist(log)) #ensure normality

meanSD <- transformed %>% 
  group_by(trait) %>% 
  summarise_at('log', .funs=list(mean=mean, sd=sd)) %>% 
  ungroup()

meanSDSpecies <- transformed %>%  
  group_by(trait, species_matched) %>% 
  summarize_at('log', .funs=list(species_mean=mean, species_sd=sd, species_length=length)) %>% 
  ungroup()

cleanContinuous <- imputedRaw %>% 
  #calculate z-scores (error risk) for continuous traits 
  left_join(transformed) %>% 
  left_join(meanSD) %>% 
  left_join(meanSDSpecies) %>% 
  mutate(error_risk_overall=(log-mean)/sd) %>% 
  mutate(error_risk_species=(log-species_mean)/species_sd) %>% 
  filter(error_risk_overall<abs(4)) %>%  #drops 195 observations (0.08% of data)
  filter(error_risk_species<abs(4)) %>% #drops an additional 2079 observations (0.90% of data)
  mutate(trait2=ifelse(trait=='leaf_area', 'Leaf Area (leaf, +petiole)',
                ifelse(trait=='SLA', 'Specific Leaf Area (+petiole)', 
                ifelse(trait=='SRL', 'Specific Root Length (all root)',
                ifelse(trait=='leaf_N', 'Leaf N Content',
                ifelse(trait=='plant_height_vegetative', 'Plant Vegetative Height',
                ifelse(trait=='seed_dry_mass', 'Seed Dry Mass',
                ifelse(trait=='leaf_dry_mass', 'Leaf Dry Mass',
                ifelse(trait=='LDMC', 'Leaf Dry Matter Content',
                trait)))))))))

# cleanContinousWide <- cleanContinuous %>% 
#   pivot_longer(cols=c('imputed_value'), names_to='data_type', values_to='trait_value') %>% 
#   na.omit()

sppNum <- cleanContinuous %>% 
  select(species_matched) %>% 
  unique()


#### Comparing original data ####

originalData <- read.csv('nutnet_trait database_combo_continuous_20240710.csv') %>% 
  select(-X3115, -X3109, -X3117, -X3114, -X614) %>% 
  pivot_longer(cols=c(plant_height_vegetative:SRL), names_to='trait', values_to='trait_value') %>% 
  filter(!is.na(trait_value)) %>%
  mutate(trait2=ifelse(trait=='leaf_area', 'Leaf Area (leaf, +petiole)',
                ifelse(trait=='SLA', 'Specific Leaf Area (+petiole)', 
                ifelse(trait=='SRL', 'Specific Root Length (all root)',
                ifelse(trait=='leaf_N', 'Leaf N Content',
                ifelse(trait=='plant_height_vegetative', 'Plant Vegetative Height',
                ifelse(trait=='seed_dry_mass', 'Seed Dry Mass',
                ifelse(trait=='leaf_dry_mass', 'Leaf Dry Mass',
                ifelse(trait=='LDMC', 'Leaf Dry Matter Content',
                trait))))))))) %>% 
  mutate(data_type=DatabaseID)

combinedContinuous <- cleanContinuous %>% 
  mutate(data_type='imputed') %>% 
  dplyr::rename(trait_value=imputed_value) %>% 
  select(DatabaseID, DatasetID, ObservationID, Family, genus, species_matched, data_type, trait_value, trait, trait2) %>% 
  rbind(originalData)



#### Boxplots for each trait ####
combinedContinuous$trait2 = factor(combinedContinuous$trait2, levels=c('Leaf Area (leaf, +petiole)', 'Leaf Dry Mass', 'Leaf Dry Matter Content', 'Specific Leaf Area (+petiole)', 'Leaf N Content', 'Plant Vegetative Height', 'Specific Root Length (all root)', 'Seed Dry Mass'))



#Look at boxplots for each trait -- means by species
combinedContinuousMean <- combinedContinuous %>%
  group_by(DatabaseID, data_type, species_matched, trait, trait2) %>% 
  dplyr::summarise(trait_value_mean=mean(trait_value)) %>% 
  ungroup()

#logged
ggplot(data=combinedContinuousMean, aes(x=as.factor(data_type), y=trait_value_mean)) +
  geom_jitter(aes(color=data_type)) +
  geom_boxplot(color='black', alpha=0) +
  facet_wrap(~trait2, scales='free_y', ncol=3, labeller=label_wrap_gen(width=25)) +
  scale_x_discrete(breaks=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY", "imputed"),
                   limits=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY", "imputed"),
                   labels=c("Au", "BN", "C2", "TP", "TY", "imp")) +
  scale_color_manual(values=c('#4E3686', '#5DA4D9', '#80D87F', 'darkgrey', '#FED23F', '#EE724C')) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position='none',
        strip.text.x = element_text(size = 20),
        axis.title.x=element_text(size=22, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=22),
        axis.title.y=element_text(size=22, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=22)) +
  xlab('Data Type') + ylab(expression(log[10]("Trait Value")))  +
  scale_y_continuous(trans='log10', labels=label_comma())
# ggsave('xyz.png', width=14, height=15, units='in', dpi=300, bg='white')

#not logged
ggplot(data=combinedContinuousMean, aes(x=as.factor(data_type), y=trait_value_mean)) +
  geom_jitter(aes(color=data_type)) +
  geom_boxplot(color='black', alpha=0) +
  facet_wrap(~trait2, scales='free_y', ncol=3, labeller=label_wrap_gen(width=25)) +
  scale_x_discrete(breaks=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY", "imputed"),
                   limits=c("AusTraits", "BIEN", "CPTD2", "TIPleaf", "TRY", "imputed"),
                   labels=c("Au", "BN", "C2", "TP", "TY", "imp")) +
  scale_color_manual(values=c('#4E3686', '#5DA4D9', '#80D87F', 'darkgrey', '#FED23F', '#EE724C')) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position='none',
        strip.text.x = element_text(size = 20),
        axis.title.x=element_text(size=22, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=22),
        axis.title.y=element_text(size=22, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=22)) +
  xlab('Data Type') + ylab("Trait Value")
# ggsave('xyz.png', width=14, height=15, units='in', dpi=300, bg='white')



##### Mean values for each species #####
meanCleanContinuous <- cleanContinuous %>% 
  group_by(Family, genus, species_matched, trait) %>% 
  dplyr::summarize(trait_value=mean(imputed_value), imputation_error=mean(std)) %>% 
  ungroup()

meanSD <- meanCleanContinuous %>% 
  mutate(log=log10(trait_value)) %>% 
  group_by(trait) %>% 
  dplyr::summarize(across('log', .fns=list(mean=mean, sd=sd))) %>% 
  ungroup()

meanSDFamily <- meanCleanContinuous %>% 
  mutate(log=log10(trait_value)) %>% 
  group_by(trait, Family) %>% 
  dplyr::summarize(across('log', .fns=list(family_mean=mean, family_sd=sd, family_length=length))) %>% 
  ungroup()

meanSDGenus <- meanCleanContinuous %>% 
  mutate(log=log10(trait_value)) %>% 
  group_by(trait, genus) %>% 
  dplyr::summarize(across('log', .fns=list(genus_mean=mean, genus_sd=sd, genus_length=length))) %>% 
  ungroup()

meanCleanContinuousErrorRisk <- meanCleanContinuous %>% 
  left_join(meanSD) %>% 
  left_join(meanSDFamily) %>% 
  left_join(meanSDGenus) %>% 
  mutate(log=log10(trait_value)) %>% 
  mutate(error_risk_overall=(log-log_mean)/log_sd, 
         error_risk_family=ifelse(log_family_length>2, (log-log_family_mean)/log_family_sd, NA),
         error_risk_genus=ifelse(log_genus_length>2, (log-log_genus_mean)/log_genus_sd, NA)) %>% 
  select(Family, genus, species_matched, trait, trait_value, imputation_error, error_risk_overall, error_risk_family, error_risk_genus) %>% 
  # left_join(meanContinuous) %>% 
  # select(-imputed_value_mean, imputed_value_sd, original_value_sd) %>% 
  mutate(trait2=ifelse(trait=='leaf_area', 'Leaf Area (leaf, +petiole)',
                ifelse(trait=='SLA', 'Specific Leaf Area (+petiole)', 
                ifelse(trait=='SRL', 'Specific Root Length (all root)',
                ifelse(trait=='leaf_N', 'Leaf N Content',
                ifelse(trait=='plant_height_vegetative', 'Plant Vegetative Height',
                ifelse(trait=='seed_dry_mass', 'Seed Dry Mass',
                ifelse(trait=='leaf_dry_mass', 'Leaf Dry Mass',
                ifelse(trait=='LDMC', 'Leaf Dry Matter Content',
                trait)))))))))


coverage <- meanCleanContinuousErrorRisk %>% 
  select(Family, species_matched, trait, trait_value) %>% 
  pivot_wider(names_from=trait, values_from=trait_value)

summary(coverage)
# 2 NAs across entire dataframe = 0.02% of data missing for these 1147 species


##### Final data #####

longContinuous <- meanCleanContinuousErrorRisk %>%
  mutate(source='Imputed Value') %>% 
  select(Family, species_matched, trait, trait_value, imputation_error, error_risk_overall, error_risk_family, error_risk_genus, source) %>% 
  dplyr::rename(species=species_matched)

# write.csv(longContinuous, 'NutNet_continuousTraitData_imputed_20240711.csv', row.names=F)
