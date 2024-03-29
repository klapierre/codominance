################################################################################
##  CoRRE_codominance.R: Calculate codominance of a community for CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: June 12, 2020
################################################################################

library(psych)
library(codyn)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\CoRRE')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
sppNames <- read.csv('corre2trykey_2021.csv')%>%
  select(genus_species, species_matched)%>%
  unique()

corre <- read.csv('CoRRE_RawAbundance_2021.csv')%>%
  # select(-X)%>%
  left_join(sppNames)%>%
  rename(old_name=genus_species, cover=abundance)%>%
  mutate(genus_species=ifelse(is.na(species_matched), as.character(old_name), as.character(species_matched)))%>%
  mutate(exp_unit=paste(site_code, project_name, community_type, plot_id, treatment, calendar_year, sep='::'))%>%
  select(-species_matched, -old_name)


##### Calculate Cmax (codominance metric) #####

#calculate relative abundance
relCover <- corre%>%
  group_by(exp_unit)%>%
  summarise(totcov=sum(cover))%>%
  ungroup()%>%
  right_join(corre)%>%
  mutate(relcov=(cover/totcov)*100)%>%
  select(-cover, -totcov)

evenness <- relCover%>%
  community_structure(time.var = 'calendar_year', abundance.var = 'relcov',
                      replicate.var = 'exp_unit', metric = c("Evar", "SimpsonEvenness", "EQ"))

# write.csv(evenness, 'corre_richEven_20240208.csv', row.names=F)

#generate rank of each species in each plot by relative cover, with rank 1 being most abundant
rankOrder <- relCover%>%
  group_by(exp_unit)%>%
  mutate(rank = as.numeric(order(order(relcov, decreasing=TRUE))))%>%
  ungroup()

###calculating harmonic means for all subsets of rank orders
#make a new dataframe with just the label
expUnit=corre%>%
  select(exp_unit)%>%
  unique()

#makes an empty dataframe
harmonicMean=data.frame(row.names=1) 

### NOTE: this code takes about 2 hours to run, so use the output in the dropbox unless there is a reason to re-run it
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
  separate(exp_unit2, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::')%>%
  mutate(calendar_year=as.integer(calendar_year))

codomSppList <- Cmax%>%
  left_join(rankOrder)%>%
  group_by(exp_unit)%>%
  filter(rank<=num_codominants)%>%
  ungroup()

# write.csv(codomSppList, 'corre_codominants_list_202402091.csv', row.names=F)

siteProjComm <- codomSppList%>%
  select(site_code, project_name, community_type)%>%
  unique()


##### Plots -- gut check if number of codominants is correct #####

rankCodominance <- Cmax %>% 
  select(exp_unit, num_codominants) %>% 
  left_join(rankOrder) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='_'))

# write.csv(rankCodominance, 'corre_codominantsRankAll_202402091.csv', row.names=F)

site_proj_comm_vector <- unique(rankCodominance$site_proj_comm)

for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(data=filter(rankCodominance, site_proj_comm == site_proj_comm_vector[PROJ]),
         aes(x=rank, y=relcov)) +
    facet_grid(rows=vars(plot_id), cols=vars(calendar_year), scales='free') +
    geom_point() +
    geom_line() +
    geom_vline(data=filter(rankCodominance, site_proj_comm == site_proj_comm_vector[PROJ]), 
               mapping=aes(xintercept=num_codominants+0.5), color="blue") +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  
  ggsave(filename=paste0("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\rank abundance curves\\",
                         site_proj_comm_vector[PROJ], "_RAC.png"),
         width = 35, height = 35, dpi = 300, units = "in", device='png')
  
}




