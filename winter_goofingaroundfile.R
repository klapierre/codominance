#Jordan goofing around file

library(tidyverse)


# Control plots. For each year, bin plots into mono-dominated, co-dominated, tri-dominated, even (4+). 
# Check to see if there is directional change over time in the control plots. 
# If yes, then decide what yr to use. If no, then use the mean (or mode?). 


#CORRE

df <- read.csv("corre_codominantsRankAll_202402091.csv")

head(df)


#playing around

#subset all data for control

controlonly <- subset(df, treatment == "C")

#ANG subset test

ANGdf <- subset(controlonly, site_code == "ANG")

unique(ANGdf$treatment_year)

#make bins

bintest <- ANGdf %>% 
  mutate(groupings = case_when(
    num_codominants == 1 ~ "monodominated",
    num_codominants == 2 ~ "codominated",
    num_codominants == 3 ~ "tridominated",
    num_codominants >= 4 ~ "even",
    TRUE ~ NA_character_
  ))

#okay i think that works how I wanted. so how to determine if it changes over time?
#unique "groupings" by plot_id?

#just eyeballing a few plots a few do change over time 
#so what year to pick??????


#do something like group by plot_id, make new column thats the avg of num_codom column, then name by like 0 -> 1.49 = monodom, 1.5 -> 2.49 = codom.....
#then do another column where we could do the same as above but do the mode. idk falling out of love with the idea of mode...


#### NEED THINK










#bin across all control

controlbinned <- controlonly %>% 
  mutate(groupings = case_when(
    num_codominants == 1 ~ "monodominated",
    num_codominants == 2 ~ "codominated",
    num_codominants == 3 ~ "tridominated",
    num_codominants >= 4 ~ "even",
    TRUE ~ NA_character_
  ))

#works, that's neat  











