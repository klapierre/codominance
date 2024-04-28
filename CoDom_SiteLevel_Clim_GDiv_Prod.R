setwd("C:/Users/alyou/OneDrive - UNCG/SIDE PROJECTS")

### LOAD PACKAGES ###
library(tidyverse)
library(codyn)


#### For 04/01/2024 ####
## CoRRE ####
CoRRERichProd <- read.csv("C:/Users/alyou/OneDrive - UNCG/2024_codominance/data/CoRRE/CoRRE_siteBiotic_2021.csv") %>%
  # keep community types w/in projects and sites separate so there are distinct ANPP and richness values for each, but the same coordinates, MAP, and MAT #
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  rename("GDiv" = "rrich", "ANPP"="anpp")
#'rrich' is gamma div#
  
CoRRECoordClim <- read.csv("C:/Users/alyou/OneDrive - UNCG/2024_codominance/data/CoRRE/CoRRE_siteLocationClimate_2021.csv") %>%
  # remove Sil and SORBAS sites #
  filter(site_code!="Sil" & site_code!="SORBAS")

CoRREfull <- merge(CoRRECoordClim,CoRRERichProd, by=c("site_code"), all=T) %>%
  select(site_code, site_proj_comm, Latitude, Longitude, MAP, MAT, GDiv, ANPP) %>%
  # remove NutNet sites from CoRRE database #
  filter(site_proj_comm!="CDR_NutNet_0" & 
           site_proj_comm!="cbgb.us_NutNet_0" &
           site_proj_comm!="shps.us_NutNet_0" &
           site_proj_comm!="sier.us_NutNet_0" &
           site_proj_comm!="temple.us_NutNet_0" &
           site_proj_comm!="veluwe.nl_NutNet_0" &
           site_proj_comm!="yarra.au_NutNet_0")

#################

## GEx ####
GExfull <- read.csv("C:/Users/alyou/OneDrive - UNCG/2024_codominance/data/GEx/GEx-metadata-with-other-env-layers-v2.csv") %>%
  rename("site_code"="site", "Latitude"="Final.Lat", "Longitude"="Final.Long", "N_deposition"="N.deposition1993",
         "MAP"="precip", "GDiv"="sprich") %>%
  # 'bio1' is (MAT*100) from WorldClim in C #
  mutate(MAT = bio1/10) %>%
  select(site_code, Latitude, Longitude, MAP, MAT, GDiv, ANPP, N_deposition) 
  # keep 'NA' for Kruger (Mananga, Shitbotawna, and Satara) sites as MAP not good predictor #
  # 'precip' is MAP #
  # 'sprich' is gamma div #

################

## NutNet ####
NutNetCoordClim <- read.csv("C:/Users/alyou/OneDrive - UNCG/2024_codominance/data/NutNet/comb-by-plot-clim-soil_2023-11-07.csv") %>%
  filter(trt=="Control") %>%
  mutate(N_Dep = as.numeric(N_Dep)) %>%
  # change site_code for 'yarra.au' since the same site_code exists in CoRRE so need to make sure they are not the same experiment #
  mutate(site_code = ifelse(site_code=="yarra.au", "yarra.au_NutNetdf", site_code)) %>%
  group_by(site_code)%>%
  summarise(Latitude=mean(latitude), Longitude=mean(longitude), GDiv=mean(site_richness),MAP=mean(MAP_v2), MAT=mean(MAT_v2), N_deposition=mean(N_Dep))
# 'site_richness' is gamma div #

NutNetprod <- read.csv("C:/Users/alyou/OneDrive - UNCG/2024_codominance/data/NutNet/full-biomass_2023-11-07.csv") %>% 
  filter(trt=="Control") %>%
  select(-year_trt, -trt, -live, -subplot, -site_name, -block) %>%
  # change site_code for 'yarra.au' since the same site_code exists in CoRRE so need to make sure they are not the same experiment (do the same for others after confirm #
  mutate(site_code = ifelse(site_code=="yarra.au", "yarra.au_NutNetdf", site_code)) %>%
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


NutNetfULL <- merge.data.frame(NutNetCoordClim, NutNetprod, by=c("site_code"), all=T)

FinalAllNetworks <- bind_rows(CoRREfull, GExfull, NutNetfULL)
# write dataframe to csv to then input to ArcGIS Pro and join with HDI raster and NDep raster data from Oak Ridge and export table to .xlsx #
write.csv(FinalAllNetworks, "C:/Users/alyou/OneDrive - UNCG/SIDE PROJECTS/\\CodominanceAllNetworkData.csv", row.names=FALSE)


