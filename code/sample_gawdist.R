
pacman::p_load(tidyverse,
               foreach)

## fake data
n_sp <- 50
df_trait <- data.frame(species = seq_len(n_sp),
                       x = runif(n_sp),
                       y = runif(n_sp),
                       z = sample(letters[1:3],
                                  size = n_sp,
                                  replace = TRUE))

df_pool <- tibble(site = rep(paste0("s", 1:3), each = 20)) %>% 
  group_by(site) %>% 
  mutate(species = sample(1:50,
                          size = n())) %>% 
  ungroup()

v_site <- unique(df_pool$site)

foreach(i = v_site) %do% {
  
  local_pool <- df_pool %>% 
    filter(site == i) %>% 
    pull(species) %>% 
    sort()

  d <- df_trait %>% 
    filter(species %in% local_pool) %>% 
    select(-species) %>% 
    gawdis()
  
  return(d)    
}
