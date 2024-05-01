################################################################################
##  codominance_drivers.R: Determining the effect of plot size and number on codominance.
##
##  Author: Kimberly Komatsu
##  Date created: January 27, 2021
################################################################################

### edited by Ashley LaRoque

# load packages
pacman::p_load(tidyverse,
               nlme,
               lsmeans,
               performance,
               ggpubr,
               DescTools)

#set working directory
#setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data') #kim's laptop

# -----ggplot theme set-----
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


# -----homemade functions-----
###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

#mode
mode <- function(codes){
  which.max(tabulate(codes))
}

# -----read in global databases (CoRRE,  GEx, NutNet)-----



#CoRRE
corre <- read_csv('data/CoRRE/corre_codominants_list_202402091.csv')%>%
  select(-block, -genus_species, -relcov, -rank)%>%
  unique()%>%
  left_join(read_csv('data/CoRRE/corre_richEven_20240208.csv'))%>%
  left_join(read_csv('data/CoRRE/corre_plot_size.csv'))%>%
  left_join(read_csv('data/CoRRE/CoRRE_siteBiotic_2021.csv'))%>%
  left_join(read_csv('data/CoRRE/CoRRE_siteLocationClimate_2021.csv'))%>%
  left_join(read_csv('data/CoRRE/CoRRE_ExperimentInfo_2021.csv'))%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(plot_number=length(plot_id))%>%
  ungroup()%>%
  mutate(database='CoRRE')%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, experiment_length, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, rrich, anpp, Cmax, num_codominants, richness, Evar)%>%
  rename(gamma_rich=rrich) %>% 
  filter(project_name != "NutNet") #remove overlapping sites from NutNet

unique(corre$trt_type)


#GEx
gex <- read_csv('data/GEx/GEx_codominants_list_20240213.csv')%>%
  select(-genus_species, -relcov, -rank)%>%
  unique()%>%
  left_join(read_csv('data/GEx/gex_richEven_20240213.csv'))%>%
  left_join(read_csv('data/GEx/GEx-metadata-with-other-env-layers-v2.csv'))%>%
  mutate(database='GEx', project_name='NA', community_type='NA', trt_type=ifelse(trt=='G', 'control', 'herb_removal'), plot_permenant='NA', MAT=bio1/10)%>%
  rename(site_code=site, plot_id=block, calendar_year=year, treatment_year=exage, plot_id=block, treatment=trt, plot_size_m2=PlotSize, MAP=bio12, gamma_rich=sprich,anpp=ANPP)%>%
  group_by(site_code, plot_id, treatment)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(plot_number=length(plot_id))%>%
  ungroup()%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)

unique(gex$trt_type)

#NutNet
nutnetANPP <- read_csv('data/nutnet/comb-by-plot-clim-soil-diversity_2023-11-07.csv')%>%
  filter(trt=='Control')%>%
  group_by(site_code, year)%>%
  summarise(anpp=mean((vascular_live_mass+nonvascular_live_mass), na.rm=T))%>%
  ungroup()%>%
  filter(!is.nan(anpp), anpp>0)%>% #drops the data from some sites in years they didn't collect, and two years at one site that reported 0 growth
  group_by(site_code)%>%
  summarise(anpp=mean(anpp, na.rm=T))%>%
  ungroup()

nutnetSiteInfo <- read_csv('data/nutnet/comb-by-plot-clim-soil-diversity_2023-11-07.csv')%>%
  group_by(site_code)%>%
  mutate(experiment_length=max(year_trt))%>%
  ungroup()%>%
  group_by(site_code, year)%>%
  mutate(plot_number=length(plot))%>%
  ungroup()%>%
  rename(MAP=MAP_v2, MAT=MAT_v2, gamma_rich=site_richness)%>%
  left_join(nutnetANPP)%>%
  select(site_code, trt, year, MAP, MAT, gamma_rich, anpp, experiment_length, plot_number)%>%
  unique()

nutnet <- read_csv('data/nutnet/NutNet_codominants_list_plot_20240213.csv')%>%
  select(exp_unit, site_code, plot, year, year_trt, trt, Cmax, num_codominants)%>%
  unique()%>%
  left_join(read_csv('data/nutnet/nutnet_plot_richEven_20240213.csv'))%>%
  left_join(nutnetSiteInfo)%>%
  mutate(database='NutNet', project_name='NA', community_type='NA', plot_size_m2=1, plot_permenant='y')%>%
  mutate(trt_type=ifelse(year_trt<1, 'control', ifelse(trt=='Control', 'control', ifelse(trt=='Fence', 'herb_removal', ifelse(trt=='NPK+Fence', 'mult_nutrient*herb_removal', ifelse(trt=='N', 'N', ifelse(trt=='P', 'P', ifelse(trt=='K', 'K', ifelse(trt=='NP', 'N*P', 'mult_nutrient')))))))))%>%
  rename(plot_id=plot, calendar_year=year, treatment_year=year_trt, treatment=trt)%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)

unique(nutnet$trt_type)

# -----combine datasets-----

#fix problem where if communities are completely even, Cmax=0 and multiple levels are listed for the plot; this needs to be fixed in original code

individualExperiments <- rbind(corre, gex, nutnet)%>%
  mutate(num_codominants_fix=ifelse(Cmax==0, richness, num_codominants))%>%
  ungroup()%>%
  select(-num_codominants)%>%
  rename(num_codominants=num_codominants_fix)%>%
  unique()

unique(individualExperiments$trt_type) # identify what treatments are

expInfo <- individualExperiments%>%
  select(exp_unit, database, site_code, project_name, community_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, trt_type)%>%
  unique()


#-----abundance cutoffs of codominance-----

correAbund <- read_csv('data/CoRRE/corre_codominantsRankAll_202402091.csv') %>% 
  select(exp_unit, site_code, project_name, community_type, plot_id, treatment, calendar_year, num_codominants, genus_species, relcov, rank)
GExAbund <- read_csv('data/GEx/gex_codominantsRankAll_202402091.csv') %>% 
  rename(site_code=site, treatment=trt, calendar_year=year) %>% 
  mutate(plot_id=paste(block, treatment, sep='_'),
         project_name=0, community_type=0) %>% 
  select(exp_unit, site_code, project_name, plot_id, community_type, treatment, calendar_year, num_codominants, genus_species, relcov, rank)
nutnetAbund <- read_csv('data/nutnet/NutNet_codominantsRankAll_20240213.csv') %>% 
  rename(calendar_year=year, plot_id=plot, treatment=trt) %>% 
  mutate(project_name=0, community_type=0) %>% 
  select(exp_unit, site_code, project_name, community_type, plot_id, treatment, calendar_year, num_codominants, genus_species, relcov, rank)

allAbund <- rbind(correAbund, GExAbund, nutnetAbund) 

replicatesAll <- allAbund %>% select(exp_unit) %>% unique() #69,886 individual data points (plot*year combinations)
replicatesSpatial <- allAbund %>% select(site_code, project_name, community_type, plot_id) %>% unique() #12,088 individual plots
replicatesExperiment <- allAbund %>% select(site_code, project_name, community_type) %>% unique() #551 experiments


#cutoff at 20% abundance for rank 1 species: 
filter1 <- allAbund %>% 
  filter(num_codominants==1 & rank==1) %>% #39,774 data points are monodominated
  mutate(drop=ifelse(relcov<30, 'c_borderline',
              ifelse(relcov<20, 'a_drop', 'b_keep'))) %>% 
  filter(drop!='b_keep') %>% #440 (1.1%) of these have cover less than 20% for the single dominant spp; 2843 (7.1%) less than 30%
  select(exp_unit, drop) %>% 
  unique()

# #cutoff at 30% abundance for sum of all codominant species -- this filter is not strict enough, don't do this
# filter2 <- allAbund %>%
#   filter(num_codominants>1) %>%
#   group_by(exp_unit) %>%
#   filter(rank<(num_codominants+1)) %>%
#   summarise(sum_cover=sum(relcov)) %>%
#   ungroup() %>%
#   filter(sum_cover<30) %>%
#   mutate(drop=1) %>%
#   select(exp_unit, drop) %>%
#   unique()
# 
# filter12 <- rbind(filter1, filter2) %>%
#   unique() #applying both filters removes 592 data points (0.85%)
# 
# filterSum <- allAbund %>% full_join(filter12) %>%
#   mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>%
#   mutate(drop2=ifelse(is.na(drop), 'keep', 'drop'))


#cutoff at 20% abundance for mean of all codominant species
filter3 <- allAbund %>% 
  filter(num_codominants>1) %>% 
  group_by(exp_unit) %>% 
  filter(rank<(num_codominants+1)) %>% 
  summarise(mean_cover=mean(relcov)) %>% 
  ungroup() %>% 
  filter(mean_cover<20) %>% 
  mutate(drop='a_drop') %>% 
  select(exp_unit, drop) %>% 
  unique()

filterMean <- rbind(filter1, filter3) %>% 
  unique() %>%  #applying both filters removes 7232 data points (10.3%)
  pivot_wider(names_from=drop, values_from=drop) %>% 
  mutate(remove=ifelse(c_borderline=='c_borderline' & a_drop=='a_drop', 1, 0)) %>% 
  filter(is.na(remove)) %>% 
  select(-remove) %>% 
  pivot_longer(c_borderline:a_drop, names_to='name', values_to='a_drop') %>% 
  filter(!is.na(a_drop)) %>% 
  select(-name) %>% 
  full_join(allAbund) %>% 
  mutate(drop2=ifelse(is.na(a_drop), 'b_keep', a_drop)) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='_'))


## new code for codom categorical rank (Ashley's edits)
## will edit 

# group number of codominants into 4 categories 
df_edit <- filterMean %>% 
  left_join(expInfo, by = c("site_code", "exp_unit", "project_name", "community_type")) %>%  # join with experimental info to understand treatments
  mutate(group = case_when(num_codominants == 1 ~ "monodominated",
                           num_codominants == 2 ~ "codominated",
                           num_codominants == 3 ~ "tridominated",
                           num_codominants >= 4 ~ "even"), # categorical version
         num_group = case_when(num_codominants == 1 ~ 1,
                          num_codominants == 2 ~ 2,
                          num_codominants == 3 ~ 3,
                          num_codominants >= 4 ~ 4)) # numerical version 

# check that there are 4 groups 
unique(df_edit$group) 

# visualize groups
ggplot(df_edit,
       aes(group)) +
  geom_bar(aes(x = factor(group, level = c('monodominated', 'codominated', 'tridominated', 'even')))) +
  theme_minimal()


#calculate mode
df_v <- df_edit %>% 
  group_by(site_proj_comm, plot_id, treatment, trt_type, calendar_year) %>% 
  reframe(mode = Mode(num_group)) %>% # needs to be checked; generates mode per group from above
  ungroup() %>% 
  mutate(treat_type = ifelse(!is.na(trt_type), trt_type, treatment)) %>% 
  filter(treat_type %in% c("control", "Control", "G")) # just control groups


df_mode <- df_v %>%  
  group_by(site_proj_comm) %>% 
  reframe(mode_yr = Mode(mode)) 
# value should be 500-900 

### Ashley is working below

ggplot(df_v,
       aes(x = mode))+
  geom_bar()

df_v2 <- df_v %>%
  pivot_wider(id_cols = c(site_proj_comm, plot_id, treatment),
              names_from = calendar_year,
              values_from = mode, names_sort = T) # wide format to see change through time 

df_v3 <- df_v %>% 
  group_by(site_proj_comm, plot_id, treatment) %>% 
  mutate( first = dplyr::first(calendar_year),
          last = dplyr::last(calendar_year),
          exp_length = last - first) %>% 
  filter(calendar_year == last)

# when trt_type = control, does trt_type = other == same 





# check control group
unique(df_edit$trt_type)

#separate control group
df_edit2 <- df_edit %>% 
  filter(trt_type == "control") # filters treatment column by control group 

edit_vector <- unique(df_edit2$site_proj_comm)

# dont run unless needed- wait time >1hr 
for(PROJ in 1:length(edit_vector)){
  ggplot(data=filter(df_edit2, site_proj_comm == edit_vector[PROJ]),
         aes(x=rank, y=relcov)) +
    geom_rect(aes(fill=group), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    scale_fill_manual(values = alpha(c("gold", "darkorange", "hotpink", "red2"), 0.2)) +
    facet_grid(rows=vars(plot_id), cols=vars(calendar_year), scales='free') +
    geom_point() +
    geom_line() +
    geom_vline(data=filter(df_edit2, site_proj_comm == edit_vector[PROJ]), 
               mapping=aes(xintercept=num_codominants+0.5), color="blue") +
    ggtitle(edit_vector[PROJ]) +
    theme_bw()
  
  ggsave(filename=paste0("/Users/ashleylaroque/Library/Mobile Documents/com~apple~CloudDocs/Grad School/Terui Lab/Codominance/codominance/fig_codom_change_",
                         edit_vector[PROJ], ".png"),
         width = 35, height = 35, dpi = 300, units = "in", device='png')
  
}

#______________________________________________


#plots for gut check
site_proj_comm_vector <- unique(filterMean$site_proj_comm)

for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(data=filter(filterMean, site_proj_comm == site_proj_comm_vector[PROJ]),
         aes(x=rank, y=relcov)) +
    geom_rect(aes(fill=drop2), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    scale_fill_manual(values = alpha(c("red", "green", "lightblue"), 0.2)) +
    facet_grid(rows=vars(plot_id), cols=vars(calendar_year), scales='free') +
    geom_point() +
    geom_line() +
    geom_vline(data=filter(filterMean, site_proj_comm == site_proj_comm_vector[PROJ]), 
               mapping=aes(xintercept=num_codominants+0.5), color="blue") +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
 
  ggsave(filename=paste0("C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_codominance\\data\\rank abundance curves\\mean cutoff\\",
                         site_proj_comm_vector[PROJ], ".png"),
         width = 35, height = 35, dpi = 300, units = "in", device='png')
  
}

#-----site-level drivers of codominance-----
nutnetControls <- individualExperiments%>%
  filter(database=='NutNet', treatment_year==0) #just pre-treatment for NutNet, so max number of plots can be included

controlsIndExp <- individualExperiments%>%
  filter(database %in% c('CoRRE', 'GEx'), trt_type=='control')%>% #control plots only
  rbind(nutnetControls)%>%
  group_by(site_code, project_name, community_type, plot_id)%>%
  summarise(num_codominants_temporal=mean(num_codominants), Evar_temporal=mean(Evar), richness_temporal=mean(richness))%>% #mean number of codominant species in a plot over time
  ungroup()%>%
  group_by(site_code, project_name, community_type)%>%
  summarise(num_codominants_mean=mean(num_codominants_temporal), num_codominants_var=var(num_codominants_temporal), Evar_mean=mean(Evar_temporal), Evar_var=var(Evar_temporal), richness_mean=mean(richness_temporal), richness_var=mean(richness_temporal))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(codominance=ifelse(num_codominants_mean<=1, 'monodominance', ifelse(num_codominants_mean>1&num_codominants_mean<=2, '2 codominants', ifelse(num_codominants_mean>2&num_codominants_mean<=3, '3 codominants', 'even'))))%>%
  mutate(num_codominants_restricted=ifelse(num_codominants_mean>5, 5, num_codominants_mean))%>%
  mutate(codom_proportion=num_codominants_mean/richness_mean)



#-----incidence of codom in corre and nutnet-----
ggplot(data=subset(controlsIndExp, database %in% c('CoRRE', 'NutNet', 'GEx')), aes(x=codominance)) +
  geom_histogram(fill='white', color='black', stat='count') +
  xlab('Number of Codominants') + ylab('Count') +
  coord_cartesian(ylim=c(0,300)) +
  scale_x_discrete(limits=c('monodominance', '2 codominants', '3 codominants', 'even'),
                   labels=c('1', '2', '3', '4+'))
#export at 400x600


#-----model - continuous co-dominance metrics-----
summary(codomSiteDrivers <- lme(num_codominants_restricted ~ MAP + MAT + gamma_rich + anpp, #10 of the 437 sites used in this analysis had num codominants reduced from >5 to 5
                                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                                random=~1|plot_size_m2))
check_model(codomSiteDrivers)
anova(codomSiteDrivers) #significant effect of MAT and gamma_rich

#R2 values
codomSiteNull <- lme(num_codominants_restricted ~ 1, 
                     data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                     random=~1|plot_size_m2)
codomMAP <- lme(num_codominants_restricted ~ MAP,
                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                random=~1|plot_size_m2)
r2(codomMAP, codomSiteNull) #MAP: marginal R2=0.000
codomMAT <- lme(num_codominants_restricted ~ MAT,
                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                random=~1|plot_size_m2)
r2(codomMAT, codomSiteNull) #MAT: marginal R2=0.023
codomGammaRich <- lme(num_codominants_restricted ~ gamma_rich,
                      data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                      random=~1|plot_size_m2)
r2(codomGammaRich, codomSiteNull) #gamma_rich: marginal R2=0.037
codomANPP <- lme(num_codominants_restricted ~ anpp,
                 data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                 random=~1|plot_size_m2)
r2(codomANPP, codomSiteNull) #ANPP: marginal R2=0.002


#-----figures of site-level drivers-----
MAPfig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                 aes(x=MAP, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('MAP (mm)') + ylab('Number of Codominants')

MATfig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                 aes(x=MAT, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('MAT (C)') + ylab('') +
  geom_smooth(method='lm', se=F, color='black', size=2)

richnessFig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                      aes(x=gamma_rich, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('Gamma Diversity') + ylab('Number of Codominants') +
  geom_smooth(method='lm', se=F, color='black', size=2)

anppFig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                  aes(x=MAP, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('Site Productivity') + ylab('')

ggarrange(MAPfig, MATfig, richnessFig, anppFig,
          ncol = 2, nrow = 2)
#export at 1200x800


#-----parameter space filled by corre and nutnet-----
ggplot(data=subset(controlsIndExp, database %in% c('CoRRE', 'NutNet')), aes(x=num_codominants_restricted)) +
  geom_density(binwidth=1, fill='white', color='black') +
  ylab('Count') + xlab('Number of Codominants')



#-----plot-level drivers of co-dominance-----
controlsPlot <- individualExperiments%>%
  mutate(keep=ifelse(database=='NutNet'&treatment_year==0, 1, ifelse(database %in% c('CoRRE', 'GEx') & trt_type=='control', 1, 0)))%>%
  filter(keep==1)%>%
  # group_by(site_code, project_name, community_type, plot_id, plot_size_m2)%>%
  # summarise(num_codominants=mean(num_codominants), richness=mean(richness), Evar=mean(Evar))%>%
  # ungroup()%>%
  mutate(num_codominants_restricted=ifelse(num_codominants>15, 15, num_codominants))%>%
  mutate(codom_proportion=num_codominants/richness)

#restricted number of codominants
summary(codomPlotDrivers <- lme(num_codominants_restricted ~ richness + Evar, 
                                data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130), 
                                random=~1|plot_size_m2))
check_model(codomPlotDrivers)
anova(codomPlotDrivers) #significant effects of plot richness and evenness

#R2 values
codomPlotNull <- lme(num_codominants_restricted ~ 1, 
                     data=subset(controlsPlot, !is.na(plot_size_m2)), 
                     random=~1|plot_size_m2)
codomRichness <- lme(num_codominants_restricted ~ richness, 
                     data=subset(controlsPlot, !is.na(plot_size_m2)), 
                     random=~1|plot_size_m2)
r2(codomRichness, codomPlotNull) #richness: marginal R2=0.042
codomEvar <- lme(num_codominants_restricted ~ Evar, 
                 data=subset(controlsPlot, !is.na(plot_size_m2)), 
                 random=~1|plot_size_m2)
r2(codomEvar, codomPlotNull) #Evar: marginal R2=0.013

#-----figures of plot-level drivers-----
richnessFig <- ggplot(data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130),
                      aes(x=richness, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('Plot Richness') + ylab('Number of Codominants') +
  geom_smooth(method='lm', se=F, color='black', size=2)

# ggplot(data=barGraphStats(data=subset(controlsPlot, !is.na(plot_size_m2)), variable="richness", byFactorNames=c("num_codominants_restricted")),
#        aes(x=as.factor(num_codominants_restricted), y=mean)) +
#   geom_bar(stat='identity') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#   xlab('Number of Codominants') + ylab('Plot Richness') +
#   coord_flip()

EvarFig <- ggplot(data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130),
                  aes(x=Evar, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('Plot Evenness') + ylab('Number of Codominants') +
  geom_smooth(method='lm', se=F, color='black', size=2)

ggarrange(richnessFig, EvarFig,
          ncol = 2, nrow = 1)
#export at 1200x600

ggplot(data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130),
       aes(x=Evar, y=num_codominants_restricted, color=richness)) +
  geom_point() +
  xlab('Plot Evenness') + ylab('Number of Codominants') +
  geom_smooth(method='lm', se=F, color='black', size=2) +
  scale_colour_gradient(trans='log', low='blue', high='red')




#-----global change treatment effects on codominance-----
correNlevels <- read_csv('CoRRE\\ExperimentInformation_March2019.csv')%>%
  filter(trt_type=='N')%>%
  select(site_code, project_name, community_type, trt_type, n)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(n_levels=length(n))%>%
  ungroup()%>%
  mutate(database='CoRRE')%>%
  select(-n)%>%
  unique()

ctlCodom <- individualExperiments%>%
  filter(treatment_year>0, trt_type=='control')%>%
  group_by(database, site_code, project_name, community_type, treatment_year, plot_size_m2)%>%
  summarise(num_codominants_control=mean(num_codominants))%>%
  ungroup()

trtCodom <- individualExperiments%>%
  filter(treatment_year>0, trt_type!='control')%>%
  group_by(database, site_code, project_name, community_type, trt_type, treatment, treatment_year, plot_size_m2)%>%
  summarise(num_codominants=mean(num_codominants))%>%
  ungroup()%>%
  left_join(ctlCodom)%>%
  mutate(codom_RR=log(num_codominants/num_codominants_control))%>%
  left_join(correNlevels)%>%
  mutate(n_levels=ifelse(database=='NutNet'&trt_type=='N', 1, ifelse(trt_type!='N', NA, n_levels)))%>%
  left_join(read_csv('CoRRE\\ExperimentInformation_March2019.csv'))%>%
  #create columns for N effect
  mutate(n=ifelse(database=='NutNet'&treatment %in% c('N', 'NP', 'NK', 'NPK', 'NPK+Fence'), 10, ifelse(database=='GEx', 0, ifelse(database=='NutNet'&treatment %in% c('P', 'K', 'PK', 'Fence', 'Control'), 0, n))))%>%
  mutate(trt_type_2=ifelse(trt_type=='N'&n<10, 'N<10', ifelse(trt_type=='N'&n>=10, 'N>10', as.character(trt_type))))%>%
  select(database, site_code, project_name, community_type, trt_type, trt_type_2, treatment, treatment_year, plot_size_m2, codom_RR, n_levels, n, p, num_codominants, num_codominants_control)%>%
  #create columns for P effect
  mutate(p=ifelse(database=='NutNet'&treatment %in% c('P', 'NP', 'PK', 'NPK', 'NPK+Fence'), 10, ifelse(database=='GEx', 0, ifelse(database=='NutNet'&treatment %in% c('N', 'K', 'NK', 'Fence', 'Control'), 0, p))))%>%
  mutate(trt_type_2=ifelse(trt_type=='P'&p<10, 'P<10', ifelse(trt_type=='P'&p>=10, 'P>10', as.character(trt_type))))%>%
  select(database, site_code, project_name, community_type, trt_type, trt_type_2, treatment, treatment_year, plot_size_m2, codom_RR, n_levels, n, p, num_codominants, num_codominants_control)%>%
  #this includes all years for each site -- consider for later analyses what to do about time
    group_by(database, site_code, project_name, community_type, trt_type, treatment, plot_size_m2)%>%
  mutate(max=max(codom_RR), min=min(codom_RR))%>%
  ungroup()%>%
  mutate(keep=ifelse(abs(min)>max, 'min', 'max'))%>%
  mutate(drop=ifelse(keep=='min' & codom_RR==min, 0, ifelse(keep=='max' & codom_RR==max, 0, 1)))%>%
  filter(drop==0)%>%
  select(-drop)

subsetTrtCodom <- trtCodom%>%
  filter(!is.na(plot_size_m2) & !is.na(trt_type) & !is.na(codom_RR) &
           trt_type %in% c('drought', 'herb_removal', 'irr', 'mult_nutrient', 'N', 'N*P', 'P', 'K', 'mult_nutrient*herb_removal'))


#-----comparing N effects in CoRRE and NutNet-----
codomN <- subsetTrtCodom%>%
  filter(trt_type=='N')%>%
  left_join(correNlevels)%>%
  mutate(database_2=ifelse(database=='NutNet', 'NutNet', ifelse(database=='CoRRE'&n<10, 'CoRRE n<10', 'CoRRE n>=10')))%>%
  mutate(database_3=ifelse(database_2 %in% c('NutNet', 'CoRRE n>=10'),  'N>=10', 'N<10'))%>%
  mutate(n_levels_cat=ifelse(n_levels>1, 'yes', 'no'))


#is there a database effect?
summary(codomNModel <- lme(codom_RR ~ as.factor(database), 
                           data=codomN, 
                           random=~1|plot_size_m2))
check_model(codomNModel)
anova(codomNModel) #CoRRE significantly lower effect than NutNet
lsmeans(codomNModel, pairwise~as.factor(database), adjust="tukey")

ggplot(data=barGraphStats(data=codomN, variable="codom_RR", byFactorNames=c("database")), aes(x=database, y=mean, fill=database)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), position=position_dodge(0.9), width=0.2) +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits=c('NutNet', 'CoRRE')) +
  scale_fill_manual(values=c('#51BBB1', '#EA8B2F')) +
  coord_cartesian(ylim=c(-0.33, 0.05)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20), legend.position='none') +
  geom_text(x=1, y=-0.33, label="*", size=10)
#export at 400x600


#is the database effect due to the level of N?
summary(codomNModel <- lme(codom_RR ~ as.factor(database_2), 
                           data=codomN, 
                           random=~1|plot_size_m2))
check_model(codomNModel)
anova(codomNModel) #no difference in N effect between CoRRE and NutNet when N added is >=10 g/m2
lsmeans(codomNModel, pairwise~as.factor(database_2), adjust="tukey")

ggplot(data=barGraphStats(data=codomN, variable="codom_RR", byFactorNames=c("database_2")), aes(x=database_2, y=mean, fill=database_2)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), position=position_dodge(0.9), width=0.2) +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits=c('NutNet', 'CoRRE n>=10', 'CoRRE n<10'),
                   labels=c('NutNet\nN=10 gm2', 'CoRRE\nN>10 gm2', 'CoRRE\nN<10 gm2')) +
  coord_cartesian(ylim=c(-0.5, 0.3)) +
  scale_fill_manual(values=c('#51BBB1', '#51BBB1', '#EA8B2F')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20), legend.position='none') +
  geom_text(x=1, y=-0.33, label="a", size=6) +
  geom_text(x=1, y=-0.34, label="   *", size=8) +
  geom_text(x=2, y=-0.45, label="a", size=6) +
  geom_text(x=2, y=-0.46, label="   *", size=8) +
  geom_text(x=3, y=0.3, label="b", size=6)
#export 600x650

#>10 gm2 vs less
ggplot(data=barGraphStats(data=codomN, variable="codom_RR", byFactorNames=c("database_3")), aes(x=database_3, y=mean, fill=database_3)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), position=position_dodge(0.9), width=0.2) +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits=c('N>=10', 'N<10'),
                   labels=c('N>10 gm2', 'N<10 gm2')) +
  coord_cartesian(ylim=c(-0.5, 0.3)) +
  scale_fill_manual(values=c('#51BBB1', 'darkgrey')) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20), legend.position='none') +
  geom_text(x=1, y=-0.34, label="   *", size=8)
#export 400x600


#N levels for threshold
summary(codomNthresholdModel <- lme(codom_RR ~ log(n), 
                                    data=subset(codomN, n_levels>1), 
                                    random=~1|plot_size_m2))
anova(codomNthresholdModel) #N level affects loss of codom spp

# # log model is best fit
# model.linear <- lm(codom_RR~n, data=subset(codomN, n_levels>1))
# model.squared <- lm(codom_RR~poly(n,2), data=subset(codomN, n_levels>1))
# model.log <- lm(codom_RR~log(n), data=subset(codomN, n_levels>1))
# anova(model.linear,model.squared)
# anova(model.linear,model.log)

# ggplot(data=subset(codomN, n_levels>1), aes(x=n, y=codom_RR, color=site_code)) +
#   geom_point() + geom_smooth(se=F, method='lm')

# ggplot(data=subset(codomN, n_levels>1), aes(x=n, y=codom_RR)) +
#   geom_point(color='#769370', size=2) +
#   geom_smooth(se=F, method='lm', color='#769370', size=2, formula='y ~ log(x)') +
#   xlab(expression(paste('Nitrogen (g',  ~ m ^ -2, ')'))) + ylab("ln RR (Number of Codominants)") +
#   geom_hline(yintercept=0) +
#   scale_x_continuous(breaks=seq(0,50,2))
# #export at 600x400

ggplot(data=subset(codomN, n_levels>1), aes(x=n, y=codom_RR)) +
  geom_point(color='darkgrey', size=3) +
  geom_smooth(se=F, method='gam', color='darkgrey', size=2) +
  xlab(expression(paste('Nitrogen (g',  ~ m ^ -2, ')'))) + ylab("ln RR (Number of Codominants)") +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-1.2,1.2))
#export at 800x600

nLevelsYes <- codomN%>%filter(n_levels>1)%>%mutate(n_levels_cat='all')
nLevelsGraph <- codomN%>%mutate(n_levels_cat=ifelse(n_levels_cat=='no', 'all', n_levels_cat))%>%rbind(nLevelsYes)

ggplot(data=subset(nLevelsGraph, database=='CoRRE'&n<60), aes(x=n, y=codom_RR, color=n_levels_cat)) +
  geom_point(size=3) +
  geom_smooth(se=F, method='gam', size=2) +
  scale_color_manual(values=c('#51BBB1','darkgrey')) +
  xlab(expression(paste('Nitrogen (g',  ~ m ^ -2, ')'))) + ylab("ln RR (Number of Codominants)") +
  geom_hline(yintercept=0) +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(-1.2,1.2))
#export at 800x600


# #random draws to illustrate gain in power from including both databases
# #make a new dataframe with just N>=10 g/m2
# codomNdrawsBoth <- codomN%>%
#   filter(database_2 %in% c('NutNet', 'CoRRE n>=10'))%>%
#   mutate(exp_unit=paste(database, site_code, project_name, community_type, treatment, sep='::'))
# 
# #makes an empty dataframe
# randomDrawsBoth=data.frame(row.names=1) 
# 
# #calculate effect size means
# for(i in 1:length(codomNdrawsBoth$exp_unit)) {
#   pull <- as.data.frame(replicate(n=1000, expr = mean(sample(codomNdrawsBoth$codom_RR, size=i, replace=F))))%>%
#     mutate(sample=i, database='Combined')
#   
#   colnames(pull)[1] <- 'mean_codom_RR'
#   
#   randomDrawsBoth=rbind(pull, randomDrawsBoth)
# }
# 
# #nutnet only
# codomNdrawsNutNet <- codomN%>%
#   filter(database_2 %in% c('NutNet'))%>%
#   mutate(exp_unit=paste(database, site_code, project_name, community_type, treatment, sep='::'))
# 
# #makes an empty dataframe
# randomDrawsNutNet=data.frame(row.names=1) 
# 
# #calculate effect size means
# for(i in 1:length(codomNdrawsNutNet$exp_unit)) {
#   pull <- as.data.frame(replicate(n=1000, expr = mean(sample(codomNdrawsNutNet$codom_RR, size=i, replace=F))))%>%
#     mutate(sample=i, database='NutNet')
#   
#   colnames(pull)[1] <- 'mean_codom_RR'
#   
#   randomDrawsNutNet=rbind(pull, randomDrawsNutNet)
# }
# 
# #corre only
# codomNdrawsCoRRE <- codomN%>%
#   filter(database_2 %in% c('CoRRE n>=10'))%>%
#   mutate(exp_unit=paste(database, site_code, project_name, community_type, treatment, sep='::'))
# 
# #makes an empty dataframe
# randomDrawsCoRRE=data.frame(row.names=1) 
# 
# #calculate effect size means
# for(i in 1:length(codomNdrawsCoRRE$exp_unit)) {
#   pull <- as.data.frame(replicate(n=1000, expr = mean(sample(codomNdrawsCoRRE$codom_RR, size=i, replace=F))))%>%
#     mutate(sample=i, database='CoRRE')
#   
#   colnames(pull)[1] <- 'mean_codom_RR'
#   
#   randomDrawsCoRRE=rbind(pull, randomDrawsCoRRE)
# }
# 
# 
# randomDraws <- rbind(randomDrawsBoth, randomDrawsNutNet, randomDrawsCoRRE)
# 
# #corre only
# ggplot(data=subset(randomDraws, database=='CoRRE'), aes(x=sample, y=mean_codom_RR, color=database)) +
#   geom_point() +
#   scale_color_manual(values=c('#51BBB1')) +
#   xlab('Sample Size') + ylab('ln RR (Number of Codominants)') +
#   geom_hline(yintercept=-0.242552, color='#51BBB1', size=2) +
#   geom_hline(yintercept=0, color='black', size=1) +
#   coord_cartesian(xlim=c(0,150))
# #export 800x600
# 
# #corre and nutnet
# ggplot(data=subset(randomDraws, database %in% c('CoRRE', 'NutNet')), aes(x=sample, y=mean_codom_RR, color=database)) +
#   geom_point() +
#   scale_color_manual(values=c('#51BBB1', '#EA8B2F')) +
#   xlab('Sample Size') + ylab('ln RR (Number of Codominants)') +
#   geom_hline(yintercept=-0.1724667, color='#EA8B2F', size=2) +
#   geom_hline(yintercept=-0.242552, color='#51BBB1', size=2) +
#   geom_hline(yintercept=0, color='black', size=1) +
#   coord_cartesian(xlim=c(0,150))
# #export 800x600
# 
# #combined
# ggplot(data=randomDraws, aes(x=sample, y=mean_codom_RR, color=database)) +
#   geom_point() +
#   scale_color_manual(values=c('black', '#51BBB1', '#EA8B2F')) +
#   xlab('Sample Size') + ylab('ln RR (Number of Codominants)') +
#   geom_hline(yintercept=-0.242552, color='#51BBB1', size=2) +
#   geom_hline(yintercept=-0.1724667, color='#EA8B2F', size=2) +
#   geom_hline(yintercept=0, color='black', size=1)  +
#   geom_hline(yintercept=-0.1913173, color='black', size=3) +
#   coord_cartesian(xlim=c(0,150))
# #export 800x600


#-----P thresholds-----
summary(codomPthresholdModel <- lme(codom_RR ~ p, 
                                    data=subset(subsetTrtCodom, trt_type=='P'), 
                                    random=~1|plot_size_m2))
anova(codomPthresholdModel) #N level affects loss of codom spp

ggplot(data=subset(subsetTrtCodom, trt_type=='P'), aes(x=p, y=codom_RR)) + geom_point() + geom_smooth(method='lm')


#-----comparing herbivore removal effects in GEx and NutNet-----
summary(codomHerbModel <- lme(codom_RR ~ as.factor(database), 
                         data=subset(subsetTrtCodom, trt_type=='herb_removal' & database %in% c('NutNet', 'GEx')), 
                         random=~1|site_code))
check_model(codomHerbModel)
anova(codomHerbModel) #no difference in herb effect between CoRRE and NutNet
lsmeans(codomHerbModel, pairwise~as.factor(database), adjust="tukey")

ggplot(data=barGraphStats(data=subset(subsetTrtCodom, trt_type=='herb_removal' & database %in% c('NutNet', 'GEx')), variable="codom_RR", byFactorNames=c("database")), aes(x=database, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#F17236') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), position=position_dodge(0.9), width=0.2) +
  scale_x_discrete(limits=c('NutNet', 'GEx')) + ylab("ln RR (Number of Codominants)") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20)) +
  geom_text(x=2, y=-0.27, label="*", size=10) 
#export at 400x600

#herbivore type
codomHerb <- subset(subsetTrtCodom, trt_type=='herb_removal' & database %in% c('NutNet', 'GEx'))%>%
  left_join(read_csv('nutnet//nutnet_grazer types.csv'))%>%
  mutate(large_grazers=ifelse(database=='GEx', 'yes', as.character(large_grazers)))

summary(codomHerbModel <- lme(codom_RR ~ as.factor(database), 
                         data=subset(codomHerb, large_grazers=='yes'), 
                         random=~1|site_code))
check_model(codomHerbModel)
anova(codomHerbModel) #no difference in herb effect between CoRRE and NutNet
lsmeans(codomHerbModel, pairwise~as.factor(database), adjust="tukey")

ggplot(data=barGraphStats(data=subset(codomHerb, large_grazers=='yes'), variable="codom_RR", byFactorNames=c("database")), aes(x=database, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#F17236') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), position=position_dodge(0.9), width=0.2) +
  scale_x_discrete(limits=c('NutNet', 'GEx')) + ylab("ln RR (Number of Codominants)") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20)) +
  geom_text(x=2, y=-0.27, label="*", size=10) 
#export at 400x600

#power together vs separate databases
with(subset(subsetTrtCodom, trt_type=='herb_removal'), t.test(codom_RR, mu=0)) #t = -2.8419, df = 280, p-value = 0.004814 ***
with(subset(subsetTrtCodom, trt_type=='herb_removal'&database=='NutNet'), t.test(codom_RR, mu=0)) # = -1.2386, df = 79, p-value = 0.2192 ***
with(subset(subsetTrtCodom, trt_type=='herb_removal'&database=='GEx'), t.test(codom_RR, mu=0)) #t = -2.9257, df = 192, p-value = 0.003851 ***


#compare across treatment types
allGCDs <- subsetTrtCodom%>%
  mutate(drop=ifelse(n<10 & trt_type %in% c('N', 'N*P', 'mult_nutrient'), 1, 0))%>% #drops any trt where N added is <10 g/m2 based on above analyses (for comparability across databases)
  filter(drop==0)%>%
  select(-drop)

summary(codomGCD <- lme(codom_RR ~ as.factor(trt_type), 
                        data=allGCDs, 
                        random=~1|site_code/project_name/community_type))
check_model(codomGCD)
anova(codomGCD) #year * trt type interaction
lsmeans(codomGCD, pairwise~as.factor(trt_type), adjust="tukey")

#difference from 0 for each treatment type
with(data=allGCDs, t.test(codom_RR, mu=0)) #t = -4.7436, df = 1019, p-value = 2.398e-06 ***
# with(subset(allGCDs, trt_type=='CO2'), t.test(codom_RR, mu=0)) #t = 2.5845, df = 10, p-value = 0.02721 ***
with(subset(allGCDs, trt_type=='drought'), t.test(codom_RR, mu=0)) #t = -0.19633, df = 22, p-value = 0.8462
with(subset(allGCDs, trt_type=='irr'), t.test(codom_RR, mu=0)) #t = -0.32187, df = 27, p-value = 0.75
# with(subset(allGCDs, trt_type=='temp'), t.test(codom_RR, mu=0)) #t = -0.88141, df = 17, p-value = 0.3904
with(subset(allGCDs, trt_type=='N'), t.test(codom_RR, mu=0)) #t = -3.5755, df = 144, p-value = 0.0004764 ***
with(subset(allGCDs, trt_type=='P'), t.test(codom_RR, mu=0)) #t = -0.79881, df = 111, p-value = 0.4261
with(subset(allGCDs, trt_type=='K'), t.test(codom_RR, mu=0)) #t = 1.0349, df = 86, p-value = 0.3036
with(subset(allGCDs, trt_type=='N*P'), t.test(codom_RR, mu=0)) #t = -2.0242, df = 118, p-value = 0.04521 ***
with(subset(allGCDs, trt_type=='mult_nutrient'), t.test(codom_RR, mu=0)) #t = -2.7949, df = 224, p-value = 0.005642 ***
with(subset(allGCDs, trt_type=='herb_removal'), t.test(codom_RR, mu=0)) #t = -2.8419, df = 280, p-value = 0.004814 ***
with(subset(allGCDs, trt_type=='mult_nutrient*herb_removal'), t.test(codom_RR, mu=0)) #t = -2.7179, df = 102, p-value = 0.007722 ***

#figure
trtBars<-allGCDs%>%
  group_by(trt_type)%>%
  summarise(mean=mean(codom_RR),
            n=length(codom_RR),
            sd=sd(codom_RR))%>%
  ungroup()%>%
  mutate(se=sd/sqrt(n))

overallBar<-allGCDs%>%
  summarise(mean=mean(codom_RR),
            n=length(codom_RR),
            sd=sd(codom_RR))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt_type='overall')

barGraph <- rbind(overallBar, trtBars)

ggplot(data=barGraph, aes(x=trt_type, y=mean, fill=trt_type)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits = c('overall', 'drought', 'irr', 'N', 'P', 'K', 'N*P', 'mult_nutrient', 'herb_removal'),
                   labels = c('overall\n(1019)', 'drt\n(22)', 'irr\n(27)', 'N\n(144)', 'P\n(111)', 'K\n(86)', 'NP\n(118)', 'mult\nnut\n(224)', 'herb\nrem\n(280)')) +
  scale_fill_manual(values=c('#F1C646', '#F17236', '#F1C646', '#769370', '#769370', '#769370', '#769370', 'black', '#769370')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  coord_cartesian(ylim=c(-0.35,0.25)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20)) +
  geom_vline(xintercept = 1.5, size = 1) +
  geom_text(x=1, y=-0.17, label="*", size=10) +
  geom_text(x=4, y=-0.32, label="*", size=10) +
  geom_text(x=7, y=-0.3, label="*", size=10) +
  geom_text(x=8, y=-0.26, label="*", size=10) +
  geom_text(x=9, y=-0.23, label="*", size=10)
#export at 1000x600




#reduced niche dimensionality
trtBars<-allGCDs%>%
  group_by(trt_type)%>%
  summarise(mean=mean(codom_RR),
            n=length(codom_RR),
            sd=sd(codom_RR))%>%
  ungroup()%>%
  mutate(se=sd/sqrt(n))

overallBar<-allGCDs%>%
  filter(trt_type %in% c('N', 'herb_removal', 'N*P', 'mult_nutrient', 'mult_nutrient*herb_removal'))%>%
  summarise(mean=mean(codom_RR),
            n=length(codom_RR),
            sd=sd(codom_RR))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt_type='overall')

barGraph <- rbind(overallBar, trtBars)

ggplot(data=subset(barGraph, trt_type %in% c('overall', 'herb_removal', 'N', 'N*P', 'mult_nutrient', 'mult_nutrient*herb_removal')), aes(x=trt_type, y=mean, fill=trt_type)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits = c('overall', 'herb_removal', 'N', 'N*P', 'mult_nutrient', 'mult_nutrient*herb_removal'),
                   labels = c('overall\n(873)', 'herb\nrem\n(280)', 'N\n(144)', 'NP\n(118)', 'mult\nnut\n(224)', 'mult nut*\nherb rem\n(103)')) +
  scale_fill_manual(values=c('#F1C646', '#769370', '#F17236', '#769370', '#769370', 'darkgrey')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  coord_cartesian(ylim=c(-0.35,0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20)) +
  geom_vline(xintercept = 1.5, size = 1) +
  geom_text(x=1, y=-0.23, label="*", size=10) +
  geom_text(x=2, y=-0.24, label="*", size=10) +
  geom_text(x=3, y=-0.33, label="*", size=10) +
  geom_text(x=4, y=-0.31, label="*", size=10) +
  geom_text(x=5, y=-0.27, label="*", size=10) +
  geom_text(x=6, y=-0.37, label="*", size=10)
#export at 1000x600


#nut added and herb removal only
ggplot(data=subset(barGraph, trt_type %in% c('herb_removal', 'N')), aes(x=trt_type, y=mean, fill=trt_type)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits = c('N', 'herb_removal'),
                   labels = c('N\n(144)', 'herb\nrem\n(280)')) +
  scale_fill_manual(values=c('#F17236', '#769370')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  coord_cartesian(ylim=c(-0.35,0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=20)) +
  geom_text(x=1, y=-0.33, label="*", size=10) +
  geom_text(x=2, y=-0.24, label="*", size=10)
#export at 400x600


#difference in codom
overall <- allGCDs%>%
  select(database, site_code, project_name, community_type, trt_type, num_codominants, num_codominants_control)%>%
  gather(key='treatment', value='num_codominants', num_codominants, num_codominants_control)%>%
  mutate(trt_ctl=ifelse(treatment=='num_codominants', 'trt', 'ctl'))

overallFig <- ggplot(data=barGraphStats(data=overall, variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'treatment\n ')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

warmFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='temp'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#F1C646') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'warming\n ')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

nFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='N'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#769370') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'N\n ')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

multNutFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='mult_nutrient'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#769370') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'multiple\nnutrients')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

herbRemFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='herb_removal'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#F17236') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'herbivore\nremoval')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

ggarrange(overallFig, warmFig, nFig, multNutFig, herbRemFig,
          ncol = 5, nrow = 1)
#export at 1200x600




#-----are the species that dominate with herb removal the same as those that dominate with N addition-----
#predict, different species win under these scenarios

#find datasets with both trt independantly
nAdded <- trtCodom%>%
  filter(n>0, trt_type %in% c('N','N*P', 'mult_nutrient'))%>%
  select(database, site_code, project_name, community_type)%>%
  unique()%>%
  mutate(N_add=1)

herbRem <- trtCodom%>%
  filter(trt_type %in% c('herb_removal'))%>%
  select(database, site_code, project_name, community_type)%>%
  unique()%>%
  mutate(herb_remove=1)

bothNaddHerbRem <- herbRem%>%
  left_join(nAdded)%>%
  mutate(both_trt=N_add+herb_remove)%>%
  filter(both_trt==2)%>%
  left_join()

sameSppCorre <- read_csv('CoRRE\\corre_codominants_list_01282021.csv')%>%
  filter(site_code %in% c('KLU', 'NIN', 'TRA'))%>%
  left_join(read_csv('CoRRE\\ExperimentInformation_March2019.csv'))%>%
  select(site_code, project_name, community_type, treatment_year, plot_id, trt_type, genus_species, relcov)%>%
  left_join(trtCodom)%>%
  filter(!is.na(database))%>%
  filter(trt_type %in% c('N', 'N*P', 'mult_nutrient', 'herb_removal'))%>%
  select(site_code, project_name, community_type, treatment_year, plot_id, trt_type, num_codominants, genus_species, relcov, codom_RR)%>%
  group_by(site_code, project_name, community_type, trt_type, num_codominants, genus_species, codom_RR)%>%
  summarise(rel_cov=mean(relcov))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, trt_type, num_codominants, codom_RR)%>%
  mutate(rank = as.numeric(order(order(rel_cov, decreasing=TRUE))))%>%
  ungroup()%>%
  filter(rank==1) #should be number of codominants?
#6 of 6 cases, different species rank 1

sameSppNutNet <- read_csv('nutnet\\NutNet_codominants_list_plot_01292021.csv')%>%
  left_join(bothNaddHerbRem)%>%
  filter(both_trt==2)%>%
  mutate(database='NutNet', project_name='NA', community_type='NA')%>%
  mutate(trt_type=ifelse(year_trt<1, 'control', ifelse(trt=='Control', 'control', ifelse(trt=='Fence', 'herb_removal', ifelse(trt=='NPK+Fence', 'mult_nutrient*herb_removal', ifelse(trt=='N', 'N', ifelse(trt=='P', 'P', ifelse(trt=='K', 'K', ifelse(trt=='NP', 'N*P', 'mult_nutrient')))))))))%>%
  rename(plot_id=plot, calendar_year=year, treatment_year=year_trt, treatment=trt)%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, num_codominants, genus_species, relcov)%>%
  select(site_code, project_name, community_type, treatment_year, plot_id, trt_type, genus_species, relcov)%>%
  left_join(trtCodom)%>%
  filter(!is.na(database))%>%
  filter(trt_type %in% c('N', 'herb_removal'))%>%
  select(site_code, project_name, community_type, treatment_year, plot_id, trt_type, num_codominants, genus_species, relcov, codom_RR)%>%
  group_by(site_code, project_name, community_type, trt_type, num_codominants, genus_species, codom_RR)%>%
  summarise(rel_cov=mean(relcov))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, trt_type, num_codominants, codom_RR)%>%
  mutate(rank = as.numeric(order(order(rel_cov, decreasing=TRUE))))%>%
  ungroup()%>%
  filter(rank==1)

herbRemNutNet <- sameSppNutNet%>%
  filter(trt_type=='herb_removal')%>%
  select(site_code, genus_species)%>%
  mutate(herb_rem=1)%>%
  unique()

nAddNutNet <- sameSppNutNet%>%
  filter(trt_type=='N')%>%
  select(site_code, genus_species)%>%
  mutate(mult_nut=1)%>%
  unique()

sameSppNutNet2 <- herbRemNutNet%>%
  full_join(nAddNutNet)%>%
  mutate(same=herb_rem+mult_nut)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  group_by(site_code)%>%
  summarise(same_spp=mean(same))%>%
  ungroup()
#42 of 72 cases, different species rank 1


#different species in N added and herb remove 48 of 78 cases = 61.5% of the time, different species dominate in the end