working.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working.dir)
library(MASS)
library(lme4)
library(effects)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ez)
library(ordinal)

#load data

d <- read.csv('clean_data/dataset.csv',header = T, sep = ',')
## compute accuracy:

d$acc <- as.numeric(d$type == d$response)
d$type <- as.character(d$type)
#accuracy summary:

prop.table(table(d$subject,d$acc),1)

ggplot(d[d$subject== c("FS","Pilot_6","pilot7"),],aes(block,rt, color = subject)) +
  geom_jitter(width = 0.05) +
  geom_boxplot(alpha = 0.3,
             color = 'black', width = 0.2) +
  geom_violin(color = 'black',
              alpha = 0.2,trim = FALSE)+
  facet_wrap(~type)

ggplot(d[d$subject== c("FS","Pilot_6","pilot7"),],aes(block, fill = as.character(acc))) +
  geom_bar(position = "fill") +
  labs(fill="Accuracy",y="Proportion")+
  scale_fill_viridis(discrete = T)+
  facet_wrap(~subject)

ggplot(d,aes(subject, fill = as.character(acc))) +
  geom_bar(position = "fill") +
  labs(fill="Accuracy",y="Proportion")+
  scale_fill_viridis(discrete = T,)


ggplot(d,aes(block,rt, color = subject)) +
  geom_jitter(width = 0.05) +
  geom_boxplot(alpha = 0.3,
               color = 'black', width = 0.2) +
  geom_violin(color = 'black',
              alpha = 0.2,trim = FALSE)+
  facet_wrap(~type)

ggplot(d[d$subject== c("FS","Pilot_6","pilot7"),],aes(block,rt, color = subject)) +
  geom_jitter(width = 0.05) +
  geom_boxplot(alpha = 0.3,
               color = 'black', width = 0.2) +
  geom_violin(color = 'black',
              alpha = 0.2,trim = FALSE)+
  facet_wrap(~type)
