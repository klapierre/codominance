################################################################################
##  NutNet_codominance.R: Calculate codominance of a community for NutNet database.
##
##  Author: Kimberly Komatsu
##  Date created: June 9, 2020
################################################################################

library(psych)
library(codyn)
library(Taxonstand)
library(WorldFlora)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\\\nutnet')

###read in data
nutnet <- read.csv('full-cover_2023-11-07.csv')%>%
  rename(cover=max_cover, genus_species=Taxon)%>%
  filter(live==1, !(genus_species %in% c('GROUND', 'OTHER LITTER', 'OTHER ARISTIDA CONTORTA (DEAD)', 'OTHER SALSOLA KALI (DEAD)', 'OTHER TRIODIA BASEDOWII (DEAD)', 'OTHER ANIMAL DROPPINGS', 'OTHER ROCK', 'OTHER ANIMAL DIGGINGS', 'OTHER WOODY OVERSTORY', 'OTHER STANDING DEAD', 'OTHER ANIMAL DIGGING', 'OTHER SOIL BIOCRUST', 'OTHER WOOD', 'OTHER SHELL', 'DEER')))%>%
  mutate(trt=as.character(ifelse(year_trt<1, 'Control', as.character(trt))))


# -----calculate Cmax (codominance metric) - plot-level-----

#calculate relative abundance
relCover <- nutnet%>%
  mutate(exp_unit=paste(site_code, block, plot, trt, year, sep='::'))%>% #group by plot
  group_by(exp_unit, site_code, block, plot, trt, year)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(nutnet)%>%
  filter(cover>0)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

evennessPlot <- relCover%>%
  community_structure(time.var = 'year', abundance.var = 'relcov',
                      replicate.var = 'exp_unit', metric = c("Evar", "SimpsonEvenness", "EQ"))

# write.csv(evennessPlot, 'nutnet_plot_richEven_20240213.csv', row.names=F)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrder <- relCover%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnit=rankOrder%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMean=data.frame(row.names=1) 

### NOTE: this code takes about 30 mins to run, so use the output in the dropbox unless there is a reason to re-run it
#calculate harmonic means
for(i in 1:length(expUnit$exp_unit)) {
  
  #creates a dataset for each unique experimental unit
  subset <- rankOrder[rankOrder$exp_unit==as.character(expUnit$exp_unit[i]),]%>%
    select(exp_unit, genus_species, relcov, rank)
  
  for(j in 1:length(subset$rank)) {
    
    #creates a dataset for each series of ranks from 1 through end of the number of ranks
    subset2 <- subset[subset$rank<=j,]
    
    #calculate harmonic mean of values
    mean <- harmonic.mean(subset2$relcov)
    meanData <- data.frame(exp_unit=unique(subset2$exp_unit),
                           num_ranks=j, 
                           harmonic_mean=mean)
    
    harmonicMean=rbind(meanData, harmonicMean)
    
  }
  
}

differenceData <- harmonicMean%>%
  left_join(rankOrder)%>%
  filter(rank==num_ranks+1)%>% #only keep the next most abundant species after the number that went into the calculation of the harmonic mean
  mutate(difference=harmonic_mean-relcov) #calculates difference between harmonic mean and the relative cover of the next most abundant species

Cmax <- differenceData%>%
  group_by(exp_unit)%>%
  summarise(Cmax=max(difference))%>%
  ungroup()%>%
  left_join(differenceData)%>%
  filter(Cmax==difference)%>%
  rename(num_codominants=num_ranks)%>%
  select(exp_unit, Cmax, num_codominants)%>%
  mutate(exp_unit2=exp_unit)%>%
  separate(exp_unit2, into=c('site', 'block', 'plot', 'trt', 'year'), sep='::')%>%
  mutate(plot=as.integer(plot), year=as.integer(year), block=as.integer(block))

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

# write.csv(codomSppList, 'NutNet_codominants_list_plot_20240213.csv', row.names=F)



##### Plots -- gut check if number of codominants is correct #####

siteProjComm <- codomSppList%>%
  select(site_code, project_name, community_type)%>%
  unique()

rankCodominance <- Cmax %>% 
  select(exp_unit, num_codominants) %>% 
  left_join(rankOrder)

# write.csv(rankCodominance, 'NutNet_codominantsRankAll_20240213.csv', row.names=F)

site_proj_comm_vector <- unique(rankCodominance$site_code)

for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(data=filter(rankCodominance, site_code == site_proj_comm_vector[PROJ]),
         aes(x=rank, y=relcov)) +
    facet_grid(rows=vars(plot), cols=vars(year), scales='free') +
    geom_point() +
    geom_line() +
    geom_vline(data=filter(rankCodominance, site_code == site_proj_comm_vector[PROJ]), 
               mapping=aes(xintercept=num_codominants+0.5), color="blue") +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  
  ggsave(filename=paste0("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\rank abundance curves\\",
                         site_proj_comm_vector[PROJ], "_RAC.png"),
         width = 35, height = 35, dpi = 300, units = "in", device='png')
  
}



# -----calculate Cmax (codominance metric) - block-level-----

#calculate relative abundance
nutnetBlock <- nutnet%>%
  mutate(exp_unit=paste(site_code, block, trt, year, sep='::'))%>% #regroup by block*trt
  group_by(exp_unit, genus_species)%>%
  summarise(cover=sum(cover), length=length(genus_species))%>%
  ungroup()

relCoverBlock <- nutnetBlock%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(nutnetBlock)%>%
  filter(cover>0)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

evennessBlock <- relCoverBlock%>%
  separate(exp_unit, c('site_code', 'block', 'trt', 'year'), sep='::')%>%
  mutate(exp_unit=paste(site_code, block, trt, sep='::'))%>%
  community_structure(time.var = 'year', abundance.var = 'relcov',
                      replicate.var = 'exp_unit', metric = c("Evar", "SimpsonEvenness", "EQ"))

# write.csv(evennessBlock, 'nutnet_block_richEven_20240213.csv', row.names=F)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrderBlock <- relCoverBlock%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnitBlock=relCoverBlock%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMeanBlock=data.frame(row.names=1) 

### NOTE: this code takes about 30 mins to run, so use the output in the dropbox unless there is a reason to re-run it
#calculate harmonic means
for(i in 1:length(expUnitBlock$exp_unit)) {
  
  #creates a dataset for each unique experimental unit
  subset <- rankOrderBlock[rankOrderBlock$exp_unit==as.character(expUnitBlock$exp_unit[i]),]%>%
    select(exp_unit, genus_species, relcov, rank)
  
  for(j in 1:length(subset$rank)) {
    
    #creates a dataset for each series of ranks from 1 through end of the number of ranks
    subset2 <- subset[subset$rank<=j,]
    
    #calculate harmonic mean of values
    mean <- harmonic.mean(subset2$relcov)
    meanData <- data.frame(exp_unit=unique(subset2$exp_unit),
                           num_ranks=j, 
                           harmonic_mean=mean)
    
    harmonicMeanBlock=rbind(meanData, harmonicMeanBlock)
    
  }
  
}

differenceDataBlock <- harmonicMeanBlock%>%
  left_join(rankOrderBlock)%>%
  filter(rank==num_ranks+1)%>% #only keep the next most abundant species after the number that went into the calculation of the harmonic mean
  mutate(difference=harmonic_mean-relcov) #calculates difference between harmonic mean and the relative cover of the next most abundant species

CmaxBlock <- differenceDataBlock%>%
  group_by(exp_unit)%>%
  summarise(Cmax=max(difference))%>%
  ungroup()%>%
  left_join(differenceDataBlock)%>%
  filter(Cmax==difference)%>%
  rename(num_codominants=num_ranks)%>%
  select(exp_unit, Cmax, num_codominants)%>%
  mutate(exp_unit2=exp_unit)%>%
  separate(exp_unit2, into=c('site', 'block', 'trt', 'year'), sep='::')%>%
  mutate(year=as.integer(year), block=as.integer(block))

codomSppListBlock <- CmaxBlock%>%
  left_join(rankOrderBlock)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

# write.csv(codomSppListBlock, 'NutNet_codominants_list_block_01292021.csv', row.names=F)



# -----calculate Cmax (codominance metric) - site-level-----
nutnetSite <- nutnet%>%
  mutate(exp_unit=paste(site_code, trt, year, sep='::'))%>% #regroup by site*trt
  group_by(exp_unit, genus_species)%>%
  summarise(cover=sum(cover), length=length(genus_species))%>%
  ungroup()

#calculate relative abundance
relCoverSite <- nutnetSite%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(nutnetSite)%>%
  filter(cover>0)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

evennessSite <- relCoverSite%>%
  separate(exp_unit, c('site_code', 'trt', 'year'), sep='::')%>%
  mutate(exp_unit=paste(site_code, trt, sep='::'))%>%
  community_structure(time.var = 'year', abundance.var = 'relcov',
                      replicate.var = 'exp_unit', metric = c("Evar", "SimpsonEvenness", "EQ"))

# write.csv(evennessSite, 'nutnet_site_richEven_01292021.csv', row.names=F)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrderSite <- relCoverSite%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnitSite=relCoverSite%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMeanSite=data.frame(row.names=1) 

### NOTE: this code takes about 30 mins to run, so use the output in the dropbox unless there is a reason to re-run it
#calculate harmonic means
for(i in 1:length(expUnitSite$exp_unit)) {
  
  #creates a dataset for each unique experimental unit
  subset <- rankOrderSite[rankOrderSite$exp_unit==as.character(expUnitSite$exp_unit[i]),]%>%
    select(exp_unit, genus_species, relcov, rank)
  
  for(j in 1:length(subset$rank)) {
    
    #creates a dataset for each series of ranks from 1 through end of the number of ranks
    subset2 <- subset[subset$rank<=j,]
    
    #calculate harmonic mean of values
    mean <- harmonic.mean(subset2$relcov)
    meanData <- data.frame(exp_unit=unique(subset2$exp_unit),
                           num_ranks=j, 
                           harmonic_mean=mean)
    
    harmonicMeanSite=rbind(meanData, harmonicMeanSite)
    
  }
  
}

differenceDataSite <- harmonicMeanSite%>%
  left_join(rankOrderSite)%>%
  filter(rank==num_ranks+1)%>% #only keep the next most abundant species after the number that went into the calculation of the harmonic mean
  mutate(difference=harmonic_mean-relcov) #calculates difference between harmonic mean and the relative cover of the next most abundant species

CmaxSite <- differenceDataSite%>%
  group_by(exp_unit)%>%
  summarise(Cmax=max(difference))%>%
  ungroup()%>%
  left_join(differenceDataSite)%>%
  filter(Cmax==difference)%>%
  rename(num_codominants=num_ranks)%>%
  select(exp_unit, Cmax, num_codominants)%>%
  mutate(exp_unit2=exp_unit)%>%
  separate(exp_unit2, into=c('site', 'trt', 'year'), sep='::')%>%
  mutate(year=as.integer(year))

codomSppListSite <- CmaxSite%>%
  left_join(rankOrderSite)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

# write.csv(codomSppListSite, 'NutNet_codominants_list_site_01292021.csv', row.names=F)



#### Species overlap with CoRRE and GEx ####

WFO.file<-read.delim("Data/CompiledData/Species_lists/WFO_Backbone/classification.txt")


nutnetSp <- TPL(unique(nutnet$genus_species))

nutnetSpClean <- nutnetSp %>% 
  filter(!is.na(New.Species), New.Species!='sp.') %>% 
  mutate(database='NutNet', species_matched=paste(New.Genus, New.Species, sep=' ')) %>% 
  select(database, species_matched) %>% 
  unique() %>% 
  mutate(species_matched=str_to_sentence(species_matched))

GExSp <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\GEx_species_family_May2023.csv') %>% 
  select(database, species_matched) %>% 
  unique()

CoRREsp <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\OriginalData\\Traits\\TRY\\corre2trykey_2021.csv') %>% 
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

needsPhotopath <- nutnet %>% 
  select(genus_species, ps_path) %>% 
  unique() %>% 
  rename(Taxon=genus_species) %>% 
  filter(ps_path=='NULL') %>% 
  left_join(nutnetSp)  %>% 
  filter(!is.na(New.Species), New.Species!='sp.') %>% 
  mutate(database='NutNet', species_matched=paste(New.Genus, New.Species, sep=' ')) %>% 
  select(species_matched, ps_path, Family) %>% 
  unique() %>% 
  mutate(species_matched=str_to_sentence(species_matched)) %>% 
  mutate(alt_photopath_possible=ifelse(Family %in% c('Acanthaceae', 'Aizoaceae', 'Amaranthaceae', 'Asteraceae', 'Boraginaceae', 'Cleomaceae', 'Caryophyllaceae', 'Cyperaceae', 'Euphorbiaceae', 'Gisekiaceae', 'Hydrocharitaceae', 'Molluginaceae', 'Nyctaginaceae', 'Polygonaceae', 'Portulacaceae', 'Poaceae', 'Scrophulariaceae', 'Zygophyllaceae', 'Cactaceae', 'Crassulaceae', 'Euphorbiaceae', 'Liliaceae', 'Bromeliaceae', 'Orchidaceae'), 'possible', 'no')) %>% 
  filter(alt_photopath_possible=='possible')

