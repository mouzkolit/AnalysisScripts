#########
#Statistical Analysis 
library(stats)
library(glm2)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(glmmADMB)
library(betareg)
stats_csv <- read.csv("fs_test.csv")
head(stats_csv, 5)
# convert as factors for further analysis
stats_csv$sex_factor <- as.factor(stats_csv$Sex)
stats_csv$treatment_factor <- as.factor(stats_csv$Treatment)
stats_csv$cno_factor <- as.factor(stats_csv$cno_first)
stats_csv$group_factor <- as.factor(stats_csv$Group)
stats_csv$mice_factor <- as.factor(stats_csv$mice_id)
stats_csv$date_factor <- as.factor(stats_csv$date_treatment)
stats_csv$effect_factor <- as.factor(stats_csv$sex_effect)

plot(stats_csv$floating_time)
hist(stats_csv$floating_time)

stats_csv$floating_dec <- stats_csv$floating_time/100


stats_out <- stats_csv[stats_csv$mice_id != 4,]
#plot boxplots and point plots 
ggplot(data = stats_csv, aes(x = treatment_factor, y = floating_time))+
  geom_boxplot()+
  geom_point(aes(color = sex_factor)) +
  geom_line(aes(group = mice_factor, color = sex_factor))
  #geom_point(position=position_jitterdodge(), aes(color = treatment_factor)) 

# check the histogram 
# Run a generalized linear mixed model since data is not 
# Normally distributed
# Negative binomial distribution
# Interaction between sex and treatment
# Mixed model with random effect using the mouse identification
fitting <- glmer(floating_dec ~ sex_factor*treatment_factor + (1|mice_factor),family="binomial", data = stats_out)

summary(fitting)
