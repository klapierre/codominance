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

#make bins test

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

### Real work

#CORRE

df <- read.csv("corre_codominantsRankAll_202402091.csv")

head(df)


#subset all data for control

controlonly <- subset(df, treatment == "C")


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

#okay so how many of the sites change between dominance states over time?

sitechangeovertime <- controlbinned %>%
  arrange(site_code, plot_id, calendar_year) %>%
  group_by(site_code, plot_id) %>%
  summarise(grouping_changes = n_distinct(groupings) > 1)

#plot dat

ggplot(sitechangeovertime, aes(x = factor(grouping_changes), fill = factor(grouping_changes))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color="black") +
  labs(x = "Grouping Changes", y = "Count") +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), 
                    labels = c("No Changes", "Changes")) +
  theme_minimal()

#okay so even in the control plots most of the plots do change between states at least once over course of exp

###
# check sandbox for resulttable1 for quick example of changes over time within a plot
###

#make column that is avg_num_codominants by plot

addavgnumcodom <- controlbinned %>%
  group_by(site_code, plot_id) %>%
  mutate(avg_num_codominants = mean(num_codominants, na.rm = TRUE)) %>%
  ungroup()

# add column that makes NEW groupings based on AVG num of codominants over time

avgnumcodomgrouping <- addavgnumcodom %>% 
  mutate(avgnum_groupings = case_when(
    between(avg_num_codominants, 0, 1.49) ~ "monodominated",
    between(avg_num_codominants, 1.5, 2.49) ~ "codominated",
    between(avg_num_codominants, 2.5, 3.49) ~ "tridominated",
    between(avg_num_codominants, 3.5, 15) ~ "even",
    TRUE ~ NA_character_
  ))









































############
# Sandbox Area
################


result_table <- addavgnumcodom %>%
  group_by(site_code, treatment_year, plot_id) %>%
  summarise(groupings = list(unique(groupings))) %>%
  ungroup()


result_table1 <- addavgnumcodom %>%
  arrange(site_code, plot_id, treatment_year) %>%
  group_by(site_code, plot_id) %>%
  mutate(treatment_sequence = row_number()) %>%
  group_by(site_code, plot_id, treatment_sequence) %>%
  summarise(treatment_year = first(treatment_year),
            groupings = list(unique(groupings))) %>%
  ungroup() %>%
  select(-treatment_sequence)

