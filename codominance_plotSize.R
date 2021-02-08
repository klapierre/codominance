################################################################################
##  codominance_plotSize.R: Determining the effect of plot size and number on codominance.
##
##  Author: Kimberly Komatsu
##  Date created: January 27, 2021
################################################################################

library(lsmeans)
library(tidyverse)

#set working directory
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\codominance\\data') #kim's laptop

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


# -----read in global experimental databases (CoRRE and GEx)-----
#CoRRE
corre <- read.csv('CoRRE\\corre_codominants_list_01282021.csv')%>%
  select(-block, -genus_species, -relcov, -rank)%>%
  unique()%>%
  left_join(read.csv('CoRRE\\corre_richEven_01292021.csv'))%>%
  left_join(read.csv('CoRRE\\corre_plot_size.csv'))%>%
  left_join(read.csv('CoRRE\\siteExperimentDetails_March2019.csv'))%>%
  left_join(read.csv('CoRRE\\ExperimentInformation_March2019.csv'))%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(plot_number=length(plot_id))%>%
  ungroup()%>%
  select(exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, rrich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)%>%
  rename(gamma_rich=rrich)

plot(corre$Evar, corre$num_codominants)
plot(corre$richness, corre$num_codominants)
plot(corre$plot_size_m2, corre$num_codominants)

#GEx
gex <- read.csv('GEx\\GEx_codominants_list_06112020.csv')%>%
  select(-genus_species, -relcov, -rank)%>%
  unique()%>%
  left_join(read.csv('GEx\\gex_richEven_01292021.csv'))%>%
  left_join(read.csv('GEx\\GEx-metadata-with-other-env-layers-v2.csv'))%>%
  mutate(project_name='NA', community_type='NA', trt_type=ifelse(trt=='G', 'control', 'herb_removal'), plot_permenant='NA', gamma_rich='NA', anpp='NA')%>%
  rename(site_code=site, plot_id=block, calendar_year=year, treatment_year=exage, plot_id=block, treatment=trt, plot_size_m2=PlotSize, MAP=bio12, MAT=bio1)%>%
  group_by(site_code, plot_id, treatment)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(plot_number=length(plot_id))%>%
  ungroup()%>%
  select(exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)

plot(gex$Evar, gex$num_codominants)
plot(gex$richness, gex$num_codominants)
plot(gex$plot_size_m2, gex$num_codominants)


# -----combine corre and gex-----
individualExperiments <- rbind(corre, gex)

expInfo <- individualExperiments%>%
  select(site_code, project_name, community_type, treatment, trt_type, plot_size_m2, plot_number, MAP, MAT, gamma_rich, anpp, experiment_length)


#-----drivers of codominance in control plots-----
controlsIndExp <- individualExperiments%>%
  filter(trt_type=='control')%>% #control plots only
  group_by(site_code, project_name, community_type, plot_id, trt_type)%>%
  summarise(num_codominants_temporal=mean(num_codominants), Evar_temporal=mean(Evar), richness_temporal=mean(richness))%>% #mean number of codominant species in a plot over time
  ungroup()%>%
  group_by(site_code, project_name, community_type, trt_type)%>%
  summarise(num_codominants_mean=mean(num_codominants_temporal), num_codominants_var=var(num_codominants_temporal), Evar_mean=mean(Evar_temporal), Evar_var=var(Evar_temporal), richness_mean=mean(richness_temporal), richness_var=mean(richness_temporal))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(codominance=ifelse(num_codominants_mean<=1.5, 'monodominance', ifelse(num_codominants_mean>1.5&num_codominants_mean<=2.5, '2 codominants', ifelse(num_codominants_mean>2.5&num_codominants_mean<=3.5, '3 codominants', 'even'))))

#model - continuous codominance metric
anova(codomPlotInfoModel <- lm(num_codominants_mean ~ plot_size_m2, data=controlsIndExp)) #plot size does not affect number of codominant species
anova(codomPlotInfoModel <- lm(num_codominants_mean ~ plot_size_m2 + Evar_mean + richness_mean, data=controlsIndExp)) #plot size does not affect number of codominant species, richness and evenness are related

# anova(evarPlotInfoModel <- lm(Evar_mean ~ plot_size_m2, data=controlsIndExp)) #plot size does not affect evenness
# anova(richPlotInfoModel <- lm(richness_mean ~ plot_size_m2, data=controlsIndExp)) #plot size does affect species richness

# ggplot(data=controlsIndExp, aes(x=plot_number, y=num_codominants_mean)) +
  # geom_point() +xlab('Number of Plots') + ylab('Number of Codominant Species')
ggplot(data=controlsIndExp, aes(x=plot_size_m2, y=num_codominants_mean)) +
  geom_point() +xlab('Plot Size (m2)') + ylab('Number of Codominant Species')

#model - categorical codominance metric
anova(lm(plot_size_m2 ~ codominance, data=controlsIndExp)) #plot size does not affect number of codominant species

ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&plot_size_m2<20), aes(x=codominance, y=plot_size_m2)) +
  geom_boxplot() +
  xlab('Number of Codominant Species') + ylab(expression(paste('Plot Size (',~m^2,')'))) +
  scale_x_discrete(limits=c('monodominance', '2 codominants', '3 codominants', 'even'),
                   labels=c('monodominance', '2 codominants', '3 codominants', '4+ codominants')) +
  coord_flip()
#export at 600x800





# -----read in coordinated global experiment database (NutNet)-----
#NutNet -- plot-level
nutnetPlot <- read.csv('nutnet\\NutNet_codominants_list_plot_01292021.csv')%>%
  select(Cmax, num_codominants, block, plot, trt, year, site_code)%>%
  unique()%>%
  rename(plot_id=plot, block_id=block)%>%
  rename(plot=num_codominants, Cmax_plot=Cmax)%>%
  group_by(site_code, trt, year)%>%
  mutate(plot_number=length(plot_id))%>%
  ungroup()

#NutNet -- block-level
nutnetBlock <- read.csv('nutnet\\NutNet_codominants_list_block_01292021.csv')%>%
  select(Cmax, num_codominants, block, trt, year, site)%>%
  unique()%>%
  rename(site_code=site, block_id=block)%>%
  rename(block=num_codominants, Cmax_block=Cmax)%>%
  group_by(site_code, trt, year)%>%
  mutate(block_number=length(block_id))%>%
  ungroup()

#NutNet -- site-level
nutnetSite <- read.csv('nutnet\\NutNet_codominants_list_site_01292021.csv')%>%
  select(Cmax, num_codominants, trt, year, site)%>%
  rename(site_code=site)%>%
  unique()%>%
  rename(site=num_codominants, Cmax_site=Cmax)

nutnetSiteInfo <- read.csv('nutnet\\comb-by-plot-clim-soil-diversity-07-December-2020.csv')%>%
  group_by(site_code, year, year_trt, trt, site_name)%>%
  summarise(MAP=mean(MAP_v2), MAT=mean(MAT_v2), gamma_rich=mean(site_richness), anpp=mean(live_mass))%>%
  ungroup()

#number of codominants
nutnetCodom <- nutnetPlot%>%
  left_join(nutnetBlock)%>%
  left_join(nutnetSite)%>%
  select(-Cmax_plot, -Cmax_block, -Cmax_site)%>%
  gather(key='scale', value='num_codominants', plot, block, site)%>%
  mutate(plot_id=ifelse(scale %in% c('block', 'site'), NA, plot_id))%>%
  mutate(block_id=ifelse(scale %in% c('site'), NA, block_id))%>%
  unique()%>%
  left_join(nutnetSiteInfo)%>%
  rename(calendar_year=year, treatment=trt, treatment_year=year_trt)%>%
  mutate(trt_type=ifelse(treatment=='Fence', 'herb_removal', ifelse(treatment=='NPK+Fence', 'mult_nutrient*herb_removal', ifelse(treatment=='Control', 'control', ifelse(treatment=='N', 'N', ifelse(treatment=='P', 'P', ifelse(treatment=='NP', 'N*P', ifelse(treatment=='K', 'K', 'mult_nutrient'))))))))%>%
  mutate(plot_size_m2=1, plot_permenant='y')%>%
  group_by(site_code)%>%
  mutate(experiment_length=length(treatment_year))%>%
  ungroup()%>%
  select(site_code, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, plot_number, block_number, scale, num_codominants)%>%
  mutate(codominance=ifelse(num_codominants<=1.5, 'monodominance', ifelse(num_codominants>1.5&num_codominants<=2.5, '2 codominants', ifelse(num_codominants>2.5&num_codominants<=3.5, '3 codominants', 'even'))))

#model - continuous codominance metric
anova(codomNutNetModel <- lm(num_codominants ~ scale, data=subset(nutnetCodom, treatment_year==0 & block_number==3 & plot_number==30))) #subset down to the "correct" NutNet model, 3 block and 30 plots per site
#plot size does affect the number of codominant species
lsmeans(codomNutNetModel, pairwise~as.factor(scale), adjust="tukey")

ggplot(data=barGraphStats(data=subset(nutnetCodom, treatment_year==0 & block_number==3 & plot_number==30), variable="num_codominants", byFactorNames=c("scale")), aes(x=scale, y=mean)) +
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Scale') + ylab('Number of Codominant Species') +
  scale_x_discrete(limits = c('plot', 'block', 'site'),
                   labels = c(expression(paste('Plot\n(1',  ~ m ^ 2, ')')),
                              expression(paste('Block\n(10',  ~ m ^ 2, ')')),
                              expression(paste('Site\n(30',  ~ m ^ 2, ')')))) +
  annotate('text', x = 1, y = 2.3, label = 'a', size = 6) +
  annotate('text', x = 2, y = 2.6, label = 'ab', size = 6) +
  annotate('text', x = 3, y = 3.05, label = 'b', size = 6) +
  theme(axis.text.x=element_text(size=16, hjust=0.5, vjust=0.5, margin=margin(t=15)))
  
#export at 600x800


# #Cmax
# nutnetCmax <- nutnetPlot%>%
#   left_join(nutnetBlock)%>%
#   left_join(nutnetSite)%>%
#   select(-plot, -block, -site)%>%
#   gather(key='scale', value='Cmax', Cmax_plot, Cmax_block, Cmax_site)%>%
#   mutate(plot_id=ifelse(scale %in% c('Cmax_block', 'Cmax_site'), NA, plot_id))%>%
#   mutate(block_id=ifelse(scale %in% c('Cmax_site'), NA, block_id))%>%
#   unique()%>%
#   left_join(nutnetSiteInfo)%>%
#   rename(calendar_year=year, treatment=trt, treatment_year=year_trt)%>%
#   mutate(trt_type=ifelse(treatment=='Fence', 'herb_removal', ifelse(treatment=='NPK+Fence', 'mult_nutrient*herb_removal', ifelse(treatment=='Control', 'control', ifelse(treatment=='N', 'N', ifelse(treatment=='P', 'P', ifelse(treatment=='NP', 'N*P', ifelse(treatment=='K', 'K', 'mult_nutrient'))))))))%>%
#   mutate(plot_size_m2=1, plot_permenant='y')%>%
#   group_by(site_code)%>%
#   mutate(experiment_length=length(treatment_year))%>%
#   ungroup()%>%
#   select(site_code, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, plot_number, block_number, scale, Cmax)
# 
# #model - continuous codominance metric
# CmaxNutNetModel <- lm(Cmax ~ scale, data=subset(nutnetCmax, treatment_year==0 & block_number==3 & plot_number==30))
# anova(CmaxNutNetModel) #plot size does not affect number of codominant species
# 
# ggplot(data=subset(nutnetCmax, treatment_year==0 & block_number==3 & plot_number>20), aes(x=scale, y=Cmax)) +
#   geom_boxplot() + 
#   xlab('Scale') + ylab('Cmax')
