################################################################################
##  database_map.R: Map of experiment locations for all databases.
##
##  Author: Kimberly Komatsu
##  Date created: February 14, 2021
################################################################################

library(ggplot2)
library(wesanderson)
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
library(tidyverse)

#set working directory
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\GEx working groups\\SEV 2019\\codominance\\data') #kim's laptop

#-----reading in data-----
corre <- read.csv('CoRRE\\siteList_LatLong.csv')%>%
  mutate(database='CoRRE (107 exp at 54 sites)')%>%
  select(-id)%>%
  rename(site_code='name')%>%
  unique()

gex <- read.csv('GEx\\GEx-metadata-with-other-env-layers-v2.csv')%>%
  mutate(database='GEx (252 sites)')%>%
  rename(site_code='site', latitude='Final.Lat', longitude='Final.Long')%>%
  select(site_code, latitude, longitude, database)%>%
  unique()
  
nutnet <- read.csv('NutNet\\comb-by-plot-clim-soil-diversity-07-December-2020.csv')%>%
  mutate(database='NutNet (117 sites)')%>%
  select(site_code, latitude, longitude, database)%>%
  unique()

latLong <- rbind(nutnet, gex, corre)


#-----making map - all three databases-----
world <- ne_countries(scale = "medium", returnclass = "sf")

z.pal <- wes_palette("Darjeeling1", 3, type = "discrete")

ggplot(data=world) +
  theme(panel.background=element_rect(fill="aliceblue", color="aliceblue")) +
  theme(text=element_text(size=20, colour="black"),
        axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black")) +
  geom_sf(color="black", fill="antiquewhite") +
  geom_point(data=latLong, mapping=aes(x=longitude, y=latitude, fill=database), size=3, shape=21) +
  scale_fill_manual(values = z.pal) +
  theme(legend.position = "top", legend.title=element_blank()) +
  ylab(element_blank()) +
  xlab(element_blank())
#export at 1000x700


#-----making map - NutNet and CoRRE only-----
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data=world) +
  theme(panel.background=element_rect(fill="white", color="white")) +
  theme(text=element_text(size=20, colour="black"),
        axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black")) +
  geom_sf(color="white", fill="lightgrey") +
  geom_point(data=subset(latLong, database %in% c('CoRRE (107 exp at 54 sites)', 'NutNet (117 sites)')), mapping=aes(x=longitude, y=latitude, fill=database), size=3, shape=21) +
  scale_fill_manual(values = c('#51BBB1', '#EA8B2F')) +
  theme(legend.position = "top", legend.title=element_blank()) +
  ylab(element_blank()) +
  xlab(element_blank())
#export at 1000x700


#-----making map - CoRRE threshold sites only-----
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data=world) +
  theme(panel.background=element_rect(fill="white", color="white")) +
  theme(text=element_text(size=20, colour="black"),
        axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black")) +
  geom_sf(color="white", fill="lightgrey") +
  geom_point(data=subset(latLong, site_code %in% c('AZI', 'CUL', 'IMGERS', 'KUFS', 'NWT', 'SVA', 'YMN')), mapping=aes(x=longitude, y=latitude), size=3, shape=21, fill='#51BBB1') +
  theme(legend.position = 'none') +
  ylab(element_blank()) +
  xlab(element_blank())
#export at 1000x700