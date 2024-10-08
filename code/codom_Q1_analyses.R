
source(here::here("code/library.R"))

# read in data frame w/ mode # of codoms #
doms <- read.csv("df_combined.csv")

  ## or for easy use in R 
df_combined <- readRDS('data_formatted/df_combined.rds')


# boxplots for each predictor #
ggplot(doms, aes(x=as.factor(mode_yr), y=abs(Latitude)))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=MAP))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=MAT))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=GDiv))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=ANPP))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=HumanDisturbance))+
  geom_boxplot()+
  coord_flip()

ggplot(doms, aes(x=as.factor(mode_yr), y=N_Deposition))+
  geom_boxplot()+
  coord_flip()

# Analyses #
a <- coef(summary(polr(as.factor(mode_yr) ~ MAP + MAT + GDiv + ANPP + HumanDisturbance + N_Deposition, data=doms, Hess=T)))
p <- pnorm(abs(a[,"t value"]), lower.tail=F)*2
b <- cbind(a, p)

