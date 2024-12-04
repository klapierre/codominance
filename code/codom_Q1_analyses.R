
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


factors <-doms[,c(10,11,12,13,14,15)]
chart.Correlation(factors, method = "spearman")

summary()

# Analyses #
a <- (multinom(as.factor(mode_yr) ~ MAP * MAT * GDiv * ANPP * HumanDisturbance * N_Deposition, data=doms))

coef <- summary(a)$coefficients
coef

stderr <- summary(a)$standard.errors
stderr

z <- coef/stderr
p_values <- 2 * (1 - pnorm(abs(z)))
p_values

