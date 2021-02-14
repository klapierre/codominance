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

nutnet <- read.csv('nutnet\\NutNet_codominants_list_plot_01292021.csv')%>%
  select(exp_unit, site_code, plot, year, year_trt, trt, Cmax, num_codominants)%>%
  unique()%>%
  left_join(read.csv('nutnet\\nutnet_plot_richEven_01292021.csv'))%>%
  left_join(nutnetSiteInfo)%>%
  mutate(database='NutNet', project_name='NA', community_type='NA', plot_size_m2=1, plot_permenant='y')%>%
  mutate(trt_type=ifelse(year_trt<1, 'control', ifelse(trt=='Control', 'control', ifelse(trt=='Fence', 'herb_removal', ifelse(trt=='NPK+Fence', 'mult_nutrient*herb_removal', ifelse(trt=='N', 'N', ifelse(trt=='P', 'P', ifelse(trt=='K', 'K', ifelse(trt=='NP', 'N*P', 'mult_nutrient')))))))))%>%
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


#-----global change treatment effects on codominance-----
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
  mutate(codom_RR=log(num_codominants/num_codominants_control))
  # mutate(drop=ifelse(database=='NutNet'&treatment %in% c('Fence', 'NPK+Fence'), 1, 0))%>% #removed NutNet fences, most of which don't really keep out herbivores in the same way that GEx exclosures do
  # filter(drop==0)%>%
  # select(-drop)

trtCodomMaxDiff <- trtCodom%>%
  group_by(database, site_code, project_name, community_type, trt_type, treatment, plot_size_m2)%>%
  mutate(max=max(codom_RR), min=min(codom_RR))%>%
  ungroup()%>%
  mutate(keep=ifelse(abs(min)>max, 'min', 'max'))%>%
  mutate(drop=ifelse(keep=='min' & codom_RR==min, 0, ifelse(keep=='max' & codom_RR==max, 0, 1)))%>%
  filter(drop==0)%>%
  select(-drop)

subsetTrtCodom <- trtCodomMaxDiff%>%
  filter(!is.na(plot_size_m2) & !is.na(trt_type) & !is.na(codom_RR) &
           trt_type %in% c('drought', 'herb_removal', 'irr', 'mult_nutrient', 'N', 'N*P', 'P', 'temp', 'K'))
  

#number of codominants RR
summary(codomGCD <- lme(codom_RR ~ as.factor(trt_type), 
                                data=subset(subsetTrtCodom, trt_type!='control'), 
                                random=~1|plot_size_m2))
check_model(codomGCD)
anova(codomGCD)
lsmeans(codomGCD, pairwise~as.factor(trt_type), adjust="tukey")

#difference from 0 for each treatment type
with(data=subsetTrtCodom, t.test(codom_RR, mu=0)) #t = -4.6149, df = 1150, p-value = 4.375e-06 ***
with(subset(subsetTrtCodom, trt_type=='drought'), t.test(codom_RR, mu=0)) #t = -0.19633, df = 22, p-value = 0.8462
with(subset(subsetTrtCodom, trt_type=='herb_removal'), t.test(codom_RR, mu=0)) #t = -2.8211, df = 280, p-value = 0.005128 ***
with(subset(subsetTrtCodom, trt_type=='irr'), t.test(codom_RR, mu=0)) #t = -0.32187, df = 27, p-value = 0.75
with(subset(subsetTrtCodom, trt_type=='mult_nutrient'), t.test(codom_RR, mu=0)) #t = -3.1666, df = 365, p-value = 0.001672 ***
with(subset(subsetTrtCodom, trt_type=='N'), t.test(codom_RR, mu=0)) #t = -2.9081, df = 183, p-value = 0.004086 ***
with(subset(subsetTrtCodom, trt_type=='N*P'), t.test(codom_RR, mu=0)) #t = -1.6492, df = 131, p-value = 0.1015 ***
with(subset(subsetTrtCodom, trt_type=='P'), t.test(codom_RR, mu=0)) #t = -0.7349, df = 111, p-value = 0.464
with(subset(subsetTrtCodom, trt_type=='temp'), t.test(codom_RR, mu=0)) #t = -0.88141, df = 17, p-value = 0.3904
with(subset(subsetTrtCodom, trt_type=='K'), t.test(codom_RR, mu=0)) #t = 0.98755, df = 86, p-value = 0.3261

#figure
trtBars<-subsetTrtCodom%>%
  group_by(trt_type)%>%
  summarise(mean=mean(codom_RR),
            n=length(codom_RR),
            sd=sd(codom_RR))%>%
  ungroup()%>%
  mutate(se=sd/sqrt(n))

overallBar<-subsetTrtCodom%>%
  summarise(mean=mean(codom_RR),
            n=length(codom_RR),
            sd=sd(codom_RR))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt_type='overall')

barGraph <- rbind(overallBar, trtBars)

ggplot(data=barGraph, aes(x=trt_type, y=mean, fill=trt_type)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("ln RR (Number of Codominants)") +
  scale_x_discrete(limits = c('overall', 'drought', 'irr', 'temp', 'N', 'P', 'K', 'N*P', 'mult_nutrient', 'herb_removal'),
                   labels = c('overall', 'drought', 'irrigation', 'warming', 'N', 'P', 'K', 'N*P', 'mult. nutrients', 'herbivore rem.')) +
  scale_fill_manual(values=c('#F1C646', '#F17236', '#F1C646', '#769370', '#769370', '#769370', '#769370', 'black', '#769370', '#F1C646')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  geom_vline(xintercept = 1.5, size = 1) +
  geom_text(x=1, y=-0.17, label="*", size=8) +
  geom_text(x=5, y=-0.23, label="*", size=8) +
  geom_text(x=9, y=-0.21, label="*", size=8) +
  geom_text(x=10, y=-0.21, label="*", size=8)


#difference in codom
overall <- subsetTrtCodom%>%
  select(database, site_code, project_name, community_type, trt_type, num_codominants, num_codominants_control)%>%
  gather(key='treatment', value='num_codominants', num_codominants, num_codominants_control)%>%
  mutate(trt_ctl=ifelse(treatment=='num_codominants', 'trt', 'ctl'))

overallFig <- ggplot(data=barGraphStats(data=overall, variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'treatment')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1))


nFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='N'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#769370') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'N')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

multNutFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='mult_nutrient'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#769370') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'mult. nutrients')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

herbRemFig <- ggplot(data=barGraphStats(data=subset(overall, trt_type=='herb_removal'), variable="num_codominants", byFactorNames=c("trt_ctl")), aes(x=trt_ctl, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#F17236') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2) +
  xlab("") +
  ylab("Number of Codominants") +
  scale_x_discrete(limits = c('ctl', 'trt'),
                   labels = c('control', 'herbivore rem.')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggarrange(overallFig, nFig, multNutFig, herbRemFig,
          ncol = 4, nrow = 1)
#export at 1200x600


#-----comparing N effects in CoRRE and NutNet-----
correNlevels <- read.csv('CoRRE\\ExperimentInformation_March2019.csv')%>%
  filter(trt_type=='N')%>%
  select(site_code, project_name, community_type, trt_type, n)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(n_levels=length(n))%>%
  ungroup()%>%
  mutate(database='CoRRE')%>%
  select(-n)%>%
  unique()

codomN <- subsetTrtCodom%>%
  filter(trt_type=='N')%>%
  left_join(correNlevels)%>%
  mutate(n_levels=ifelse(database=='NutNet', 1, n_levels))%>%
  left_join(read.csv('CoRRE\\ExperimentInformation_March2019.csv'))%>%
  select(database, site_code, project_name, community_type, treatment_year, n_levels, n, plot_size_m2, codom_RR, num_codominants, num_codominants_control)%>%
  mutate(n=ifelse(database=='NutNet', 10, n))

#is there a database effect? - not significnat, but a trend for stronger responses in NutNet
summary(codomNModel <- lme(codom_RR ~ as.factor(database), 
                        data=codomN, 
                        random=~1|plot_size_m2))
check_model(codomNModel)
anova(codomNModel) #no difference in N effect between CoRRE and NutNet
lsmeans(codomNModel, pairwise~as.factor(database), adjust="tukey")


ggplot(data=barGraphStats(data=codomN, variable="codom_RR", byFactorNames=c("database")), aes(x=database, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill='#769370') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab("Database") +
  ylab("ln RR (Number of Codominants)")

#N amount
summary(codomNlevelModel <- lme(codom_RR ~ n, 
                      data=codomN, 
                      random=~1|plot_size_m2))
anova(codomNlevelModel) #no difference in N effect between CoRRE and NutNet

ggplot(data=codomN, aes(x=n, y=codom_RR, color=database)) +
  geom_point() + geom_smooth(se=F, method='lm')

#N levels for threshold
summary(codomNthresholdModel <- lme(codom_RR ~ n, 
                                data=subset(codomN, n_levels>1), 
                                random=~1|plot_size_m2))
anova(codomNthresholdModel) #no difference in N effect between CoRRE and NutNet

ggplot(data=subset(codomN, n_levels>1), aes(x=n, y=codom_RR, color=site_code)) +
  geom_point() + geom_smooth(se=F, method='lm', formula = y ~ x + I(x^2))


#-----comparing herbivore removal effects in GEx and NutNet-----
summary(codomHerb <- lme(codom_RR ~ as.factor(database), 
                        data=subset(subsetTrtCodom, trt_type=='herb_removal' & database %in% c('NutNet', 'GEx')), 
                        random=~1|site_code))
check_model(codomHerb)
anova(codomHerb) #no difference in N effect between CoRRE and NutNet
lsmeans(codomHerb, pairwise~as.factor(database), adjust="tukey")

ggplot(data=barGraphStats(data=subset(subsetTrtCodom, trt_type=='herb_removal' & database %in% c('NutNet', 'GEx')), variable="codom_RR", byFactorNames=c("database")), aes(x=database, y=mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  xlab("Database") +
  ylab("ln RR (Number of Codominants)")
