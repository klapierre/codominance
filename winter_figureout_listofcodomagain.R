df <- df_grouped %>% ungroup()

### Do this by year and by plot?
### anytime there are a pair/trio of species designated codom, keep for
### comparisons? 
### Pull out ANYTHING, that was paired as a codom.

### List of codominating sp, for every plot and year where a 
# codominant pairing EVER occurs

# group by site, projectid,year, plot - sort alphabetically,
# cast it wide, get unique

codom_mode <- df %>% 
  group_by(site_code, project_name, community_type, plot_id) %>% 
  summarise(plotmode = mode(num_codominants)) %>% 
  ungroup() %>% 
  group_by(site_code, project_name, community_type) %>% 
  summarise(mode = mode(plotmode)) %>% 
  ungroup() 

df_fix <- df %>% 
  left_join(codom_mode) %>% 
  #filter(mode >1, mode <4) %>% 
  filter(rank <= num_codominants) # %>% 
  # filter(num_codominants == mode)


?str_equal  
