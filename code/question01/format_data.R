# load all data from corre, nutnet, and gex
# merge together and clean

source("code/library.R")
source(here::here("code/functions.R"))

# categorical groups of codoms: details on line 206
  # df_grouped<- readRDS("data_formatted/df_grouped.rds")

# mode of codoms: details on line 233
  # df_mode_q1 <- readRDS("data_formatted/df_mode_q1.rds")


# Read data ---------------------------------------------------------------


corre <- read_csv('data/corre/1corre_codominants_list_202402091.csv')%>%
  dplyr::select(-block, -genus_species, -relcov, -rank)%>%
  unique()%>%
  left_join(read_csv('data/corre/corre_richEven_20240208.csv'))%>%
  left_join(read_csv('data/corre/corre_plot_size.csv'))%>%
  left_join(read_csv('data/corre/corre_siteBiotic_2021.csv'))%>%
  left_join(read_csv('data/corre/corre_siteLocationClimate_2021.csv'))%>%
  left_join(read_csv('data/corre/corre_ExperimentInfo_2021.csv')) %>% 
  group_by(site_code, project_name, community_type, calendar_year, treatment) %>%
  mutate(plot_number=length(plot_id),
         database='corre') %>%
  ungroup() %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length = max(treatment_year)) %>%
  ungroup() %>%
  dplyr::select(database, exp_unit, site_code, project_name, community_type, plot_id, 
                calendar_year, treatment_year, experiment_length, treatment, trt_type,
                plot_size_m2, plot_number, plot_permenant, MAP, MAT, rrich, anpp, Cmax,
                num_codominants, richness, Evar) %>%
  rename(gamma_rich=rrich) %>% 
  filter(project_name != "NutNet")


unique(corre$trt_type)


#gex

gex <- read_csv('data/gex/1gex_codominants_list_20240213.csv')%>%
  dplyr::select(-genus_species, -relcov, -rank)%>%
  unique()%>% 
  left_join(read_csv('data/gex/gex_richEven_20240213.csv'))%>%
  left_join(read_csv('data/gex/gex-metadata-with-other-env-layers-v2.csv'))%>%
  unique() %>% 
  mutate(database='gex', 
         project_name='NA', 
         community_type='NA',
         trt_type=ifelse(trt=='G', 'control', 'herb_removal'), 
         plot_permenant='NA', 
         MAT=bio1/10)%>%
  rename(site_code = site,
         plot_id = block, 
         calendar_year = year, 
         treatment_year = exage,
         plot_id = block, 
         treatment = trt, 
         plot_size_m2 = PlotSize, 
         MAP = bio12,
         gamma_rich = sprich,
         anpp = ANPP)%>%
  group_by(site_code, plot_id, treatment) %>%
  mutate(experiment_length = max(treatment_year)) %>%
  ungroup() %>%
  group_by(site_code, project_name, community_type, calendar_year, treatment) %>%
  mutate(plot_number = length(plot_id)) %>%
  ungroup() %>%
  dplyr::select(database, exp_unit, site_code, project_name, community_type, plot_id, 
                calendar_year, treatment_year, treatment, trt_type, plot_size_m2, 
                plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length,
                Cmax, num_codominants, richness, Evar)

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
  dplyr::select(site_code, trt, year, MAP, MAT, gamma_rich, anpp, experiment_length, plot_number)%>%
  unique()

nutnet <- read_csv('data/nutnet/NutNet_codominants_list_plot_20240213.csv')%>%
  dplyr::select(exp_unit, site_code, plot, year, year_trt, trt, Cmax, num_codominants)%>%
  unique()%>%
  left_join(read_csv('data/nutnet/nutnet_plot_richEven_20240213.csv')) %>% 
  left_join(nutnetSiteInfo) %>%
  mutate(database='NutNet', project_name='NA', community_type='NA', plot_size_m2=1, plot_permenant='y',
         trt_type=ifelse(year_trt<1, 'control',
                         ifelse(trt=='Control', 'control',
                                ifelse(trt=='Fence', 'herb_removal', 
                                       ifelse(trt=='NPK+Fence', 'mult_nutrient*herb_removal',
                                              ifelse(trt=='N', 'N',
                                                     ifelse(trt=='P', 'P', 
                                                            ifelse(trt=='K', 'K', 
                                                                   ifelse(trt=='NP', 'N*P', 'mult_nutrient')))))))))%>%
  rename(plot_id=plot, calendar_year=year, treatment_year=year_trt, treatment=trt)%>%
  dplyr::select(database, exp_unit, site_code, project_name, community_type, plot_id,
                calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number,
                plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, 
                num_codominants, richness, Evar)

unique(nutnet$trt_type)

# -----combine datasets-----

#fix problem where if communities are completely even, Cmax=0 and multiple levels are listed for the plot; this needs to be fixed in original code

individualExperiments <- rbind(corre, gex, nutnet) %>%
  mutate(num_codominants_fix = ifelse(Cmax==0, richness, num_codominants)) %>%
  ungroup() %>%
  dplyr::select(-num_codominants) %>%
  rename(num_codominants = num_codominants_fix) %>%
  unique()

unique(individualExperiments$trt_type) # identify what treatments are

expInfo <- individualExperiments%>%
  dplyr::select(exp_unit, database, site_code, project_name, community_type, 
                plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp,
                trt_type)%>%
  unique()


#-----abundance cutoffs of codominance-----

correAbund <- read_csv('data/corre/rank_corre_codominants_202402091.csv') %>% 
  dplyr::select(exp_unit, site_code, project_name, community_type, plot_id, treatment, calendar_year, num_codominants, genus_species, relcov, rank)

gexAbund <- read_csv('data/GEx/gex_codominantsRankAll_202402091.csv') %>% 
  rename(site_code=site, treatment=trt, calendar_year=year) %>% 
  mutate(plot_id=paste(block, treatment, sep='_'),
         project_name=0, community_type=0) %>% 
  dplyr::select(exp_unit, site_code, project_name, plot_id, community_type, treatment, calendar_year, num_codominants, genus_species, relcov, rank)

nutnetAbund <- read_csv('data/nutnet/NutNet_codominantsRankAll_20240213.csv') %>% 
  rename(calendar_year=year, plot_id=plot, treatment=trt) %>% 
  mutate(project_name=0, community_type=0) %>% 
  dplyr::select(exp_unit, site_code, project_name, community_type, plot_id, treatment, calendar_year, num_codominants, genus_species, relcov, rank)

allAbund <- rbind(correAbund, gexAbund, nutnetAbund) 

replicatesAll <- allAbund %>% dplyr::select(exp_unit) %>% unique() #69,886 individual data points (plot*year combinations)
replicatesSpatial <- allAbund %>% dplyr::select(site_code, project_name, community_type, plot_id) %>% unique() #12,088 individual plots
replicatesExperiment <- allAbund %>% dplyr::select(site_code, project_name, community_type) %>% unique() #551 experiments


#cutoff at 20% abundance for rank 1 species: 
filter1 <- allAbund %>% 
  filter(num_codominants==1 & rank==1) %>% #39,774 data points are monodominated
  mutate(drop=ifelse(relcov<30, 'c_borderline',
                     ifelse(relcov<20, 'a_drop', 'b_keep'))) %>% 
  filter(drop!='b_keep') %>% #440 (1.1%) of these have cover less than 20% for the single dominant spp; 2843 (7.1%) less than 30%
  dplyr::select(exp_unit, drop) %>% 
  unique()


#cutoff at 20% abundance for mean of all codominant species
filter3 <- allAbund %>% 
  filter(num_codominants>1) %>% 
  group_by(exp_unit) %>% 
  filter(rank<(num_codominants+1)) %>% 
  summarise(mean_cover=mean(relcov)) %>% 
  ungroup() %>% 
  filter(mean_cover<20) %>% 
  mutate(drop='a_drop') %>% 
  dplyr::select(exp_unit, drop) %>% 
  unique()

filterMean <- rbind(filter1, filter3) %>% 
  unique() %>%  #applying both filters removes 7232 data points (10.3%)
  pivot_wider(names_from=drop, values_from=drop) %>% 
  mutate(remove=ifelse(c_borderline=='c_borderline' & a_drop=='a_drop', 1, 0)) %>% 
  filter(is.na(remove)) %>% 
  dplyr::select(-remove) %>% 
  pivot_longer(c_borderline:a_drop, names_to='name', values_to='a_drop') %>% 
  filter(!is.na(a_drop)) %>% 
  dplyr::select(-name) %>% 
  full_join(allAbund) %>% 
  mutate(drop2=ifelse(is.na(a_drop), 'b_keep', a_drop)) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='_'))


# Create categorical codom groups -----------------------------------------

# group number of codominants into 4 categories 
df_grouped <- filterMean %>% 
  left_join(expInfo, by = c("site_code", "exp_unit", "project_name", "community_type")) %>%  # join with experimental info to understand treatments
  mutate(group = case_when(num_codominants == 1 ~ "monodominated",
                           num_codominants == 2 ~ "codominated",
                           num_codominants == 3 ~ "tridominated",
                           num_codominants >= 4 ~ "even"),
         num_group = case_when(num_codominants == 1 ~ 1,
                               num_codominants == 2 ~ 2,
                               num_codominants == 3 ~ 3,
                               num_codominants >= 4 ~ 4)) %>% 
  group_by(site_proj_comm) %>% 
  mutate(mean_dominants_raw = mean(num_codominants), # mean number of dominants from raw classification
         mean_dominacnts_mod = mean(num_group)) # mean number of dominants from modified groups

saveRDS(df_grouped, file = "data_formatted/df_grouped.rds") # only added columns to simplify # of codoms

# visualize groups
ggplot(df_grouped,
       aes(group)) +
  geom_bar(aes(x = factor(group, level = c('monodominated', 'codominated', 'tridominated', 'even')))) +
  theme_minimal()



# Calculate mode ----------------------------------------------------------


# calculate mode inc control for each year, site, proj, community, treatment, plot
df_mode <- df_grouped %>%  
  mutate(treat_type = ifelse(!is.na(trt_type), trt_type, treatment)) %>% 
  group_by(site_proj_comm, site_code, project_name, community_type, plot_id, treat_type) %>% 
  reframe(plot_codom = Mode(num_group)) %>% # mode function must be capital here 
  ungroup()  

# subset controls 
df_control <- df_mode %>%
  filter(treat_type %in% c("control", "Control", "G"))
# above: is the treatment "reference" also a control group?
# above: there are some other items in 'treatment' that are labeled as control in 'trt_type'/'treat_type', is this correct?

# calculate mode across all years of a treatment just for control groups 
df_mode_q1 <- df_control %>%  
  group_by(site_code, project_name, community_type) %>% # mode generated from these
  reframe(mode_yr = Mode(plot_codom)) %>%  
  ungroup()

saveRDS(df_mode_q1, file = "data_formatted/df_mode_q1.rds")


