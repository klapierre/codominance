setwd("C:/Users/alyoung6/OneDrive - UNCG/SIDE PROJECTS")

### LOAD PACKAGES ###
library(multcomp)
library(tidyverse)
library(codyn)
library(vegan)
library(abdiv)
library(lme4)
library(lmerTest)
library(emmeans)

#### For 04/01/2024 ####
## CoRRE ####
CoRRERichProd <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/CoRRE/CoRRE_siteBiotic_2021.csv")%>%
  # keep community types w/in projects and sites separate so there are distinct ANPP and richness values for each, but the same coordinates, MAP, and MAT #
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  rename("GDiv" = "rrich", "ANPP"="anpp")
#'rrich' is gamma div#

CoRRECoordClim <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/CoRRE/CoRRE_siteLocationClimate_2021.csv") %>%
  # remove Sil and SORBAS sites #
  filter(site_code!="Sil" & site_code!="SORBAS")

CoRREfull <- merge(CoRRECoordClim,CoRRERichProd, by=c("site_code"), all=T) %>%
  select(site_code, site_proj_comm, Latitude, Longitude, MAP, MAT, GDiv, ANPP)

#################

## GEx ####
GExfull <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/GEx/GEx-metadata-with-other-env-layers-v2.csv") %>%
  rename("site_code"="site", "Latitude"="Final.Lat", "Longitude"="Final.Long", "N_deposition"="N.deposition1993",
         "MAP"="precip", "GDiv"="sprich") %>%
  select(site_code, Latitude, Longitude, MAP, GDiv, ANPP, N_deposition) 
  # keep 'NA' for Kruger (Mananga, Shitbotawna, and Satara) sites as MAP not good predictor #
  # 'precip' is MAP #
  # 'sprich' is gamma div #

################

## NutNet ####
NutNetCoordClim <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/NutNet/comb-by-plot-clim-soil_2023-11-07.csv") %>%
  mutate(N_Dep = as.numeric(N_Dep)) %>%
  group_by(site_code)%>%
  summarise(Latitude=mean(latitude), Longitude=mean(longitude), GDiv=mean(site_richness),MAP=mean(MAP_v2), MAT=mean(MAT_v2), N_deposition=mean(N_Dep))
# 'site_richness' is gamma div #

NutNetprod <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/NutNet/full-biomass_2023-11-07.csv") %>% 
  filter(trt=="Control") %>%
  select(-year_trt, -trt, -live, -subplot, -site_name, -block) %>%
  group_by(year, site_code, plot) %>%
  spread(category,mass) %>%
  rename("FORB_Phlox_diffusa" = "FORB + PHLOX DIFFUSA") %>%
  group_by(year, site_code, plot) %>%
  # Ingrid said "LIVE" is all live biomass that doesn't fit into one of the other categories #
  mutate(TOTALLY = sum(GRAMINOID, WOODY, FORB, LEGUME, PTERIDOPHYTE, VASCULAR, LIVE, ANNUAL, PERENNIAL, BRYOPHYTE, CACTUS, FORB_Phlox_diffusa, LICHEN,na.rm=T)) %>%
  gather(category, mass, 4:22) %>%
  group_by(year, site_code, plot) %>%
  summarise(anpp = ifelse(category=="TOTAL"|category=="TOTALLY", mass, NA)) %>%
  # remove negative ANPP value from twostep.us before averaging #
  filter(anpp > 0) %>%
  group_by(site_code) %>%
  summarise(ANPP = mean(anpp))


NutNetfULL <- merge.data.frame(NutNetClim, NutNetproductivity, by=c("site_code"), all=T)

FinalAllNetworks <- bind_rows(CoRREfull, GExfull, NutNetfULL)
write.csv(FinalAllNetworks, "C:/Users/alyoung6/OneDrive - UNCG/SIDE PROJECTS/\\CodominanceAllNetworksSummaryAbioticRichProd.csv", row.names=FALSE)


################################################################################################################################
#### For 03/25/2024 ####
# Read in data for CoRRE #
CoRREcodom <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/CoRRE/corre_codominantsRankAll_202402091.csv")
CoRRElocat <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/CoRRE/CoRRE_siteLocationClimate_2021.csv")
CoRRE <- merge(CoRREcodom, CoRRElocat, by=c("site_code"), all=T) %>%
  rename("year"="calendar_year") %>%
  select(year, genus_species, site_code, Latitude, Longitude)
# Read in data for GEx #
GExcodom <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/GEx/gex_codominantsRankAll_202402091.csv")
GExlocat <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/GEx/GEx-metadata-with-other-env-layers-v2.csv")
GEx <- merge(GExcodom, GExlocat, by=c("site"), all=T) %>%
  rename("site_code"="site", "Latitude"="Final.Lat", "Longitude"="Final.Long") %>%
  select(year, genus_species, site_code, Latitude, Longitude)


# Read in data for NutNet #
NutNetcodom <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/NutNet/NutNet_codominantsRankAll_20240213.csv") %>%
  select(site_code, year, genus_species)
NutNetlocat <- read.csv("C:/Users/alyoung6/OneDrive - UNCG/2024_codominance/data/NutNet/comb-by-plot_2023-11-07.csv") %>%
  group_by(site_code)%>%
  summarise(Latitude=mean(latitude), Longitude=mean(longitude))%>%
  select(site_code, Latitude, Longitude)
NutNet <- merge(NutNetcodom, NutNetlocat, by=c("site_code"), all=T) %>%
  select(year, genus_species, site_code, Latitude, Longitude)


FinalAllNetworks <- rbind(CoRRE, GEx, NutNet)
write.csv(FinalAllNetworks, "C:/Users/alyoung6/OneDrive - UNCG/SIDE PROJECTS/\\CodominanceAllNetworkSppLatLon.csv", row.names=FALSE)
#############################################################################################################################################