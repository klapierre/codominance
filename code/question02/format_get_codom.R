
# setup -------------------------------------------------------------------

pacman::p_load(DescTools,
               tidyverse)

# data --------------------------------------------------------------------

## - cleaned species keys
## - need to confirm with Kim
df_key_nutnet <- read_csv("data/nutnet/NutNet_clean_spp_names_20240710.csv",
                          guess_max = 10000) %>% 
  mutate(latin_name = ifelse(is.na(New.Species), 
                             NA,
                             paste(New.Genus, New.Species))) %>% 
  select(genus_species = Taxon, 
         latin_name)

df_key_corre <- read_csv("data/corre/trykey_corre2_2021.csv") %>% 
  select(genus_species,
         latin_name = species_matched)

df_key_gex <- read_csv("data/gex/gex_sppfam_final_10june2020.csv") %>% 
  select(genus_species,
         latin_name = clean_ejf)

## - combine species keys across datasets
df_key <-  bind_rows(df_key_nutnet,
                     df_key_corre,
                     df_key_gex) %>% 
  distinct(genus_species,
           latin_name)

## - join df0
df0 <- readRDS("data_formatted/df_grouped.rds") %>%
  ungroup() %>% 
  left_join(df_key,
            by = "genus_species")

# get mode ----------------------------------------------------------------

## DescTools::Mode() can return more than one value when there are ties
## take max() after DescTools::Mode() to make it a scalar
## QUESTION: does plot_id accounts for treatments? Are there any swap in treament within a given plot?
df_mode <- df0 %>%
  group_by(site_code,
           project_name, 
           community_type, 
           plot_id) %>% 
  summarise(mode = Mode(num_codominants) %>% 
              max(),
            .groups = "drop")

## - retain only those of # codominance = 2 or 3
## - get the most common pair/trio over time
df_codom <- df0 %>% 
  mutate(cutoff = ifelse(num_codominants < 4,
                         yes = num_codominants,
                         no = 4)) %>% 
  group_by(site_code,
           project_name,
           community_type,
           plot_id,
           calendar_year) %>% 
  mutate(r = rank(-relcov)) %>%
  ungroup() %>% 
  filter(r <= cutoff) %>% 
  distinct(site_code,
           project_name,
           community_type,
           plot_id,
           genus_species) %>% 
  arrange(site_code, 
          project_name, 
          community_type, 
          plot_id) %>% 
  group_by(site_code,
           project_name,
           community_type,
           plot_id) %>% 
  mutate(n_sp = n_distinct(genus_species)) %>% 
  ungroup()

## test visualization - should be removed later
df_codom %>% 
  distinct(site_code, 
           project_name, 
           community_type,
           plot_id,
           n_sp) %>% 
  ggplot(aes(x = n_sp)) +
  geom_histogram(binwidth = 0.5) +
  labs(x = "Number of co-dominant species (at least once in the observation period)")

# #join to the starting df0
# df_fix <- df0 %>%
#   left_join(codom_mode) %>%
#   #filter(mode >1, mode <4) %>%
#   filter(rank <= num_codominants) # %>%
# # filter(num_codominants == mode)
# 
# #make the nice list (the arrange might be unnecessary but I'm ripping from old code so I'm not getting adventurous)
# nice_df <- df_fix %>%
#   select(site_code, project_name, community_type, calendar_year, plot_id,num_codominants, genus_species) %>%
#   group_by(site_code, project_name, community_type, calendar_year, plot_id) %>%
#   arrange(site_code, project_name, community_type, calendar_year, plot_id, num_codominants, genus_species) %>%
#   mutate(row_number = row_number()) %>%
#   ungroup() %>%
#   pivot_wider(names_from = genus_species, values_from = row_number, values_fn = length, values_fill = list(row_number = 0)) %>%
#   filter(num_codominants %in% c(2,3))
# 
# 
# ### ^^^ This works I think! ^^^ ####
# 
# ### Part 2 ###
# # Get a complete species list for each site_code, project_name, community_type. #
# 
# species_list <- df_fix %>%
#   group_by(site_code, project_name, community_type) %>%
#   summarise(unique_genus_species = list(unique(genus_species)), .groups = 'drop')
# 
# # okay so this one makes a row entry that is a concatenated list of all the unique. do we want this a different way? #
# # But it works though lol #
# 
# 
# #### Below are notes from our discussion: Can ignore ###
# 
# ### Do this by year and by plot?
# ### anytime there are a pair/trio of species designated codom, keep for
# ### comparisons?
# ### Pull out ANYTHING, that was paired as a codom.
# 
# ### List of codominating sp, for every plot and year where a
# # codominant pairing EVER occurs
# 
# # group by site, projectid,year, plot - sort alphabetically,
# # cast it wide, get unique
# 
# 
# 
# 
# 
# 
# ?str_equal
