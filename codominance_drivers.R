################################################################################
##  codominance_drivers.R: Determining the effect of plot size and number on codominance.
##
##  Author: Kimberly Komatsu
##  Date created: January 27, 2021
################################################################################

library(nlme)
library(lsmeans)
library(performance)
library(ggpubr)
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

#mode
mode <- function(codes){
  which.max(tabulate(codes))
}

# -----read in global databases (CoRRE,  GEx, NutNet)-----
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
  mutate(database='CoRRE')%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, rrich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)%>%
  rename(gamma_rich=rrich)

#GEx
gex <- read.csv('GEx\\GEx_codominants_list_06112020.csv')%>%
  select(-genus_species, -relcov, -rank)%>%
  unique()%>%
  left_join(read.csv('GEx\\gex_richEven_01292021.csv'))%>%
  left_join(read.csv('GEx\\GEx-metadata-with-other-env-layers-v2.csv'))%>%
  mutate(database='GEx', project_name='NA', community_type='NA', trt_type=ifelse(trt=='G', 'control', 'herb_removal'), plot_permenant='NA', MAT=bio1/10)%>%
  rename(site_code=site, plot_id=block, calendar_year=year, treatment_year=exage, plot_id=block, treatment=trt, plot_size_m2=PlotSize, MAP=bio12, gamma_rich=sprich,anpp=ANPP)%>%
  group_by(site_code, plot_id, treatment)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(plot_number=length(plot_id))%>%
  ungroup()%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)

#NutNet
nutnetANPP <- read.csv('nutnet\\comb-by-plot-clim-soil-diversity-07-December-2020.csv')%>%
  filter(trt=='Control')%>%
  group_by(site_code, year)%>%
  summarise(anpp=mean(live_mass, na.rm=T))%>%
  ungroup()%>%
  filter(!is.nan(anpp), anpp>0)%>% #drops the data from some sites in years they didn't collect, and two years at one site that reported 0 growth
  group_by(site_code)%>%
  summarise(anpp=mean(anpp, na.rm=T))%>%
  ungroup()

nutnetSiteInfo <- read.csv('nutnet\\comb-by-plot-clim-soil-diversity-07-December-2020.csv')%>%
  mutate(database='NutNet', project_name='NA', community_type='NA', plot_size_m2=1, plot_permenant='y')%>%
  mutate(trt_type=ifelse(trt=='Control', 'control', ifelse(trt=='Fence', 'herb_removal', ifelse(trt=='NPK+Fence', 'mult_nutrient*herb_removal', ifelse(trt=='N', 'N', ifelse(trt=='P', 'P', ifelse(trt=='K', 'K', ifelse(trt=='NP', 'N*P', 'mult_nutrient'))))))))%>%
  group_by(site_code)%>%
  mutate(experiment_length=max(year_trt))%>%
  ungroup()%>%
  group_by(site_code, year)%>%
  mutate(plot_number=length(plot))%>%
  ungroup()%>%
  rename(MAP=MAP_v2, MAT=MAT_v2, gamma_rich=site_richness)%>%
  left_join(nutnetANPP)%>%
  select(database, site_code, project_name, community_type, trt, trt_type, year, MAP, MAT, gamma_rich, anpp, experiment_length, plot_number, plot_size_m2, plot_permenant)%>%
  unique()

nutnet <- read.csv('nutnet\\NutNet_codominants_list_plot_01292021.csv')%>%
  select(exp_unit, site_code, plot, year, year_trt, trt, Cmax, num_codominants)%>%
  unique()%>%
  left_join(read.csv('nutnet\\nutnet_plot_richEven_01292021.csv'))%>%
  left_join(nutnetSiteInfo)%>%
  rename(plot_id=plot, calendar_year=year, treatment_year=year_trt, treatment=trt)%>%
  select(database, exp_unit, site_code, project_name, community_type, plot_id, calendar_year, treatment_year, treatment, trt_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp, experiment_length, Cmax, num_codominants, richness, Evar)


# -----combine datasets-----
individualExperiments <- rbind(corre, gex, nutnet)

expInfo <- individualExperiments%>%
  select(database, site_code, project_name, community_type, plot_size_m2, plot_number, plot_permenant, MAP, MAT, gamma_rich, anpp)%>%
  unique()


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
  mutate(codominance=ifelse(num_codominants_mean<=1.5, 'monodominance', ifelse(num_codominants_mean>1.5&num_codominants_mean<=2.5, '2 codominants', ifelse(num_codominants_mean>2.5&num_codominants_mean<=3.5, '3 codominants', 'even'))))%>%
  mutate(num_codominants_restricted=ifelse(num_codominants_mean>5, 5, num_codominants_mean))%>%
  mutate(codom_proportion=num_codominants_mean/richness_mean)



#-----model - continuous co-dominance metrics-----
summary(codomSiteDrivers <- lme(num_codominants_restricted ~ MAP + MAT + gamma_rich + anpp, #10 of the 437 sites used in this analysis had num codominants reduced from >5 to 5
                           data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                           random=~1|plot_size_m2))
check_model(codomSiteDrivers)
anova(codomSiteDrivers)

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
r2(codomMAT, codomSiteNull) #MAT: marginal R2=0.024
codomGammaRich <- lme(num_codominants_restricted ~ gamma_rich,
                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                random=~1|plot_size_m2)
r2(codomGammaRich, codomSiteNull) #gamma_rich: marginal R2=0.044
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
  xlab('MAT (C)') + ylab('Number of Codominants') +
  geom_smooth(method='lm', se=F, color='black', size=2)

richnessFig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                      aes(x=gamma_rich, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('Gamma Diversity') + ylab('Number of Codominants') +
  geom_smooth(method='lm', se=F, color='black', size=2)

anppFig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                  aes(x=MAP, y=num_codominants_restricted)) +
  geom_point(color='grey45') +
  xlab('Site Productivity') + ylab('Number of Codominants')

ggarrange(MAPfig, MATfig, richnessFig, anppFig,
          ncol = 2, nrow = 2)
#export at 1200x800


#-----model - proportional co-dominance-----
summary(codomSiteDrivers <- lme(codom_proportion ~ MAP + MAT + gamma_rich + anpp, #10 of the 437 sites used in this analysis had num codominants reduced from >5 to 5
                                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                                random=~1|plot_size_m2))
check_model(codomSiteDrivers)
anova(codomSiteDrivers)

#R2 values
codomSiteNull <- lme(codom_proportion ~ 1, 
                     data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                     random=~1|plot_size_m2)
codomMAP <- lme(codom_proportion ~ MAP,
                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                random=~1|plot_size_m2)
r2(codomMAP, codomSiteNull) #MAP: marginal R2=0.018
codomMAT <- lme(codom_proportion ~ MAT,
                data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                random=~1|plot_size_m2)
r2(codomMAT, codomSiteNull) #MAT: marginal R2=0.005
codomGammaRich <- lme(codom_proportion ~ gamma_rich,
                      data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                      random=~1|plot_size_m2)
r2(codomGammaRich, codomSiteNull) #gamma_rich: marginal R2=0.167
codomANPP <- lme(codom_proportion ~ anpp,
                 data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)), 
                 random=~1|plot_size_m2)
r2(codomANPP, codomSiteNull) #ANPP: marginal R2=0.100


#-----figures of site-level drivers-----
MAPfig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                 aes(x=MAP, y=codom_proportion)) +
  geom_point(color='grey45') +
  xlab('MAP (mm)') + ylab('Number of Codominants/Plot Richness') +
  geom_smooth(method='lm', se=F, color='black', size=2)

MATfig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                 aes(x=MAT, y=codom_proportion)) +
  geom_point(color='grey45') +
  xlab('MAT (C)') + ylab('Number of Codominants/Plot Richness') +
  geom_smooth(method='lm', se=F, color='black', size=2)

richnessFig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                      aes(x=gamma_rich, y=codom_proportion)) +
  geom_point(color='grey45') +
  xlab('Gamma Diversity') + ylab('Number of Codominants/Plot Richness') +
  geom_smooth(method='lm', se=F, color='black', size=2)

anppFig <- ggplot(data=subset(controlsIndExp, !is.na(plot_size_m2)&!is.na(MAP)&!is.na(MAT)&!is.na(gamma_rich)&!is.na(anpp)),
                  aes(x=MAP, y=codom_proportion)) +
  geom_point(color='grey45') +
  xlab('Site Productivity') + ylab('Number of Codominants/Plot Richness') +
  geom_smooth(method='lm', se=F, color='black', size=2)

ggarrange(MAPfig, MATfig, richnessFig, anppFig,
          ncol = 2, nrow = 2)
#export at 1200x800


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
anova(codomPlotDrivers)

#R2 values
codomPlotNull <- lme(num_codominants_restricted ~ 1, 
                     data=subset(controlsPlot, !is.na(plot_size_m2)), 
                     random=~1|plot_size_m2)
codomRichness <- lme(num_codominants_restricted ~ richness, 
                      data=subset(controlsPlot, !is.na(plot_size_m2)), 
                      random=~1|plot_size_m2)
r2(codomRichness, codomPlotNull) #richness: marginal R2=0.041
codomEvar <- lme(num_codominants_restricted ~ Evar, 
                 data=subset(controlsPlot, !is.na(plot_size_m2)), 
                 random=~1|plot_size_m2)
r2(codomEvar, codomPlotNull) #Evar: marginal R2=0.014

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



# #ratio of codominants to richness
# summary(codomPlotDrivers <- lme(codom_proportion ~ richness + Evar, 
#                                 data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130), 
#                                 random=~1|plot_size_m2))
# check_model(codomPlotDrivers)
# anova(codomPlotDrivers)
# 
# #R2 values
# codomPlotNull <- lme(codom_proportion ~ 1, 
#                      data=subset(controlsPlot, !is.na(plot_size_m2)), 
#                      random=~1|plot_size_m2)
# codomRichness <- lme(codom_proportion ~ richness, 
#                      data=subset(controlsPlot, !is.na(plot_size_m2)), 
#                      random=~1|plot_size_m2)
# r2(codomRichness, codomPlotNull) #richness: marginal R2=0.133
# codomEvar <- lme(codom_proportion ~ Evar, 
#                  data=subset(controlsPlot, !is.na(plot_size_m2)), 
#                  random=~1|plot_size_m2)
# r2(codomEvar, codomPlotNull) #Evar: marginal R2=0.033
# 
# 
# #-----figures of plot-level drivers-----
# richnessFig <- ggplot(data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130),
#                       aes(x=richness, y=codom_proportion)) +
#   geom_point(color='grey45') +
#   xlab('Plot Richness') + ylab('Number of Codominants') +
#   geom_smooth(method='lm', se=F, color='black', size=2)
# 
# eVarFig <- ggplot(data=subset(controlsPlot, !is.na(plot_size_m2) & richness<130),
#                   aes(x=Evar, y=codom_proportion)) +
#   geom_point(color='grey45') +
#   xlab('Plot Evenness') + ylab('Number of Codominants/Plot Richness') +
#   geom_smooth(method='lm', se=F, color='black', size=2)
# 
# 
# ggarrange(MAPfig, MATfig, richnessFig, anppFig,
#           ncol = 2, nrow = 2)
# #export at 1200x800