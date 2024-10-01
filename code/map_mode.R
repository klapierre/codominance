library(tidyverse)
library(DescTools)

df_mode_q1 <- readRDS("data_formatted/df_mode_q1.rds")

df_mode_ashley <- df_mode_q1 %>% 
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep= "_"))


Codominance_AllSiteData <- read_csv("Codominance_AllSiteData.csv")
Codominance_AllSiteData <- Codominance_AllSiteData %>% 
  mutate(site_proj_comm = ifelse(is.na(site_proj_comm), site_code, site_proj_comm))

df_combined <- df_mode_ashley %>% 
  left_join(Codominance_AllSiteData, by = "site_proj_comm")











corre_df <- read.csv("corre_codominantsRankAll_202402091.csv")
gex_codominantsRankAll_202402091 <- read_csv("gex_codominantsRankAll_202402091.csv")
NutNet_codominantsRankAll_20240213 <- read_csv("NutNet_codominantsRankAll_20240213.csv")




# Brought in site level biotic and abiotic data from CoBuddies Group


controlonlyCoRRe <- subset(corre_df, treatment %in% c("C",
                                                      "c",
                                                      "N0F0",
                                                      "u_u_c",
                                                      "Camb_Namb",
                                                      "Control",
                                                      "control",
                                                      "CNA",
                                                      "0",
                                                      "CONTROL",
                                                      "mixed_CO",
                                                      "t1",
                                                      "9",
                                                      "9_f_u_n",
                                                      "N0S0H0",
                                                      "N0P0S0",
                                                      "con",
                                                      "CA_N1",
                                                      "1_0_CO",
                                                      "A",
                                                      "uuuu",
                                                      "_000",
                                                      "CURRENT",
                                                      "amb",
                                                      "Open_Ungrazed",
                                                      "1NF",
                                                      "N0B0",
                                                      "ref_rich16",
                                                      "gcc",
                                                      "N0",
                                                      "CT",
                                                      "AcAt",
                                                      "N0M0",
                                                      "x",
                                                      "N0P0",
                                                      "0_CONTROL",
                                                      "CK_W0",
                                                      "ambient",
                                                      "ct",
                                                      "CK",
                                                      "N1P0",
                                                      "ambient_control",
                                                      "Cont",
                                                      "XXX",
                                                      "P3N0",
                                                      "T0F0",
                                                      "OO",
                                                      "UnwarmedControl",
                                                      "Reference",
                                                      "CC",
                                                      "U CC",
                                                      "0_CONTROL_1",
                                                      "0N0P")) %>% 
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>% 
  # remove Sil and SORBAS sites #
  filter(site_code!="Sil" & site_code!="SORBAS") %>% 
  # remove NutNet sites from CoRRE database #
  filter(site_proj_comm!="CDR_NutNet_0" & 
           site_proj_comm!="cbgb.us_NutNet_0" &
           site_proj_comm!="shps.us_NutNet_0" &
           site_proj_comm!="sier.us_NutNet_0" &
           site_proj_comm!="temple.us_NutNet_0" &
           site_proj_comm!="veluwe.nl_NutNet_0" &
           site_proj_comm!="yarra.au_NutNet_0")

controlonlyCoRRe$block <- as.character(controlonlyCoRRe$block)


controlonlyGex <- subset(gex_codominantsRankAll_202402091, trt %in% "U") %>% 
  group_by(site, year) %>% 
  #summarise(AVG_Site_Mode = mean(num_codominants)) %>% 
  mutate(site_proj_comm = site) %>% 
  ungroup()#averaging across blocks within a site within a year


controlnlyNut <- subset(NutNet_codominantsRankAll_20240213, trt %in% "Control") %>% 
  mutate(site_proj_comm = site_code)
controlnlyNut$block <- as.character(controlnlyNut$block)

MegaBind <- bind_rows(controlnlyNut,controlonlyCoRRe,controlonlyGex)  %>% 
  #mutate(site_proj_comm = ifelse(site_proj_comm == NA, site_code, site_proj_comm)) 
  mutate(site_proj_comm = ifelse(is.na(site_proj_comm), site, site_proj_comm))  


ModeCodom <- MegaBind %>%
  group_by(site_proj_comm) %>% 
  summarise(SiteDomMode = Mode(num_codominants)) %>% 
  rename(site_proj_comm = "site_proj_comm")



MapData <- left_join(Codominance_AllSiteData, ModeCodom, by = "site_proj_comm")

write.csv(MapData,"MapData.csv")

NewData <- MapData %>% 
  mutate(Doms = ifelse(SiteDomMode == 1, "Monodominated", ifelse(SiteDomMode == 2, "Codominated", "Even"))) %>% 
  na.omit(MapData$SiteDomMode)
NewData$Doms <- factor(NewData$Doms, levels = c("Monodominated", "Codominated", "Even"))

autumnalPalette <- c("#02385A", "#A63922", "#D8B573")

NewData %>% 
  ggplot(aes(x="", y=SiteDomMode, fill = Doms)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y",start=0)

ggplot(NewData, aes(x = Doms, fill = Doms)) +
  geom_histogram(stat = "count") +
  stat_count(binwidth = 1, 
             geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5))+
  xlab("Site Dominance")+
  ylab("Count")+
  scale_fill_manual(values = autumnalPalette) +
  theme(legend.position = "none")
