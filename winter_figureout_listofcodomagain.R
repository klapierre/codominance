library(tidyverse)

### Part 1 ###
# Get a nice list of the codominant species for each plot and year that has 2 or 3 codominants ###

#load this from the data_formatted folder in the github
df <- df_grouped %>% ungroup()

#get the codom mode 
codom_mode <- df %>% 
  group_by(site_code, project_name, community_type, plot_id) %>% 
  summarise(plotmode = mode(num_codominants)) %>% 
  ungroup() %>% 
  group_by(site_code, project_name, community_type) %>% 
  summarise(mode = mode(plotmode)) %>% 
  ungroup() 

#join to the starting df
df_fix <- df %>% 
  left_join(codom_mode) %>% 
  #filter(mode >1, mode <4) %>% 
  filter(rank <= num_codominants) # %>% 
# filter(num_codominants == mode)

#make the nice list (the arrange might be unnecessary but I'm ripping from old code so I'm not getting adventurous)
nice_df <- df_fix %>%
  select(site_code, project_name, community_type, calendar_year, plot_id,num_codominants, genus_species) %>% 
  group_by(site_code, project_name, community_type, calendar_year, plot_id) %>% 
  arrange(site_code, project_name, community_type, calendar_year, plot_id, num_codominants, genus_species) %>% 
  mutate(row_number = row_number()) %>%  
  ungroup() %>%
  pivot_wider(names_from = genus_species, values_from = row_number, values_fn = length, values_fill = list(row_number = 0)) %>% 
  filter(num_codominants %in% c(2,3))


### ^^^ This works I think! ^^^ ####
 
### Part 2 ###
# Get a complete species list for each site_code, project_name, community_type. #

species_list <- df_fix %>% 
  group_by(site_code, project_name, community_type) %>%
  summarise(unique_genus_species = list(unique(genus_species)), .groups = 'drop')

# okay so this one makes a row entry that is a concatenated list of all the unique. do we want this a different way? #
# But it works though lol #


#### Below are notes from our discussion: Can ignore ###
 
### Do this by year and by plot?
### anytime there are a pair/trio of species designated codom, keep for
### comparisons? 
### Pull out ANYTHING, that was paired as a codom.

### List of codominating sp, for every plot and year where a 
# codominant pairing EVER occurs

# group by site, projectid,year, plot - sort alphabetically,
# cast it wide, get unique






?str_equal  
