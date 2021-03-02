
# #install.packages('ggplot2')
# install.packages('reshape')
# install.packages('lme4')
# install.packages('ez')
# install.packages('emmeans')
# install.packages('gridExtra')
# install.packages('knitr')

library(ggplot2)
library(reshape)
library(lme4)
library(ez)
library(emmeans)
library(gridExtra)
library(knitr)


# filename <- '/Users/iris.mencke/Documents/Peak_MGAs_four_channels.csv'

filename <- '/Volumes/Projects/2017-0121-CCMusic/analysis_MM/Peak_MGAs_four_channels.csv'

d <- read.csv(filename, header = TRUE, sep = ",", dec = ".",quote = "'")
d <- subset(d, select = -c(atonal_comps,tonal_comps) )
d_long <- melt(d, id.vars = 1:2, measure.vars = seq(3,34,2))
d_lat <- melt(d, id.vars = 1:2, measure.vars = seq(4,34,2))
d_long$lat <- d_lat$value
xx <- t(data.frame(strsplit(as.character(d_long$variable), "_")))
d_long <- cbind(d_long, xx[,c(1,2,4)])

rownames(d_long) <- c()
colnames(d_long)[c(4,6,7,8)] <- c('amp', 'hem', 'dev', 'cond')
d_long$variable <- c() # delete column

#d_long$amp[grepl('left', d_long$hem)] <- d_long$amp[grepl('left', d_long$hem)]*-1

d_long$cond <- factor(d_long$cond, levels = c('tonal', 'atonal'))
d_long$dev <- factor(d_long$dev, levels = c('pitch', 'location', 'intensity', 'timbre'))

########################################################################
## AMPLITUDES
#######################################################################


# 1. fit mixed effect models

uamp0 <- lmer(amp~1+(1|ID), data = d_long, REML = FALSE) # model without any effects; null-model
uamp1 <- lmer(amp~cond+(1|ID), data = d_long, REML = FALSE) # model with effect of condition; + random effect; random intercept change the basline per subject to mimic variance in a population
uamp2 <- lmer(amp~cond+dev + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
uamp3 <- lmer(amp~cond+dev+hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
uamp4 <- lmer(amp~cond+dev+hem + cond:dev + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
uamp5 <- lmer(amp~cond+dev+hem + cond:dev + cond:hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
uamp6 <- lmer(amp~cond+dev+hem + cond*dev + cond*hem + dev*hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
uamp7 <- lmer(amp~cond+dev+hem + cond:dev + cond:hem + dev:hem + cond:dev:hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond

# Likelihood ratio test and store models:

uamps.test <- anova(uamp0,uamp1,uamp2,uamp3,uamp4,uamp5,uamp6,uamp7); 
uamps.test

uamps.models <- list(uamp0,uamp1,uamp2,uamp3,uamp4,uamp5,uamp6,uamp7);

# 2. PAIRWISE CONTRASTS

# TABLE 2: pairwise comparisons with three factors

confint.uamps3 <- as.data.frame(confint(lsmeans(uamp7,pairwise~cond|dev|hem,
                                                adjust = "Bonferroni")$contrast))

pairwise.uamps3 <-cbind(as.data.frame(lsmeans(uamp7, pairwise~cond|dev|hem,
                                              adjust="bonferroni")$contrasts),
                        confint.uamps3[c("lower.CL","upper.CL")])

pairwise.uamps3 <- pairwise.uamps3[c("contrast","dev","hem","estimate","lower.CL",
                                     "upper.CL","t.ratio","p.value")]

pairwise.uamps3[,"Cohen's d"] <- pairwise.uamps3$estimate/sigma(uamp7)
pairwise.uamps3[,4:8] <- round(pairwise.uamps3[,4:8],2)
colnames(pairwise.uamps3) <- c("contrast","deviant","hemisphere","estimate",
                               "CI 2.5%","CI 97.5%","t","p","d")
kable(pairwise.uamps3)

# make table
write.table(pairwise.uamps3,file= "/Volumes/Projects/2017-0121-CCMusic/analysis_MM/pairwise_cond_uamp7.csv",
            sep=",",row.names = FALSE,quote=FALSE)

###################################################################################

# TABLE 4a: Pairwise contrasts of mean amplitudes between features.

confint.uamps.hem <- as.data.frame(confint(lsmeans(uamp7,pairwise~dev,
                                                   adjust = "Bonferroni")$contrast))

pairwise.uamps.hem <-cbind(as.data.frame(lsmeans(uamp7, pairwise~dev,
                                                 adjust="bonferroni")$contrasts),
                           confint.uamps.hem[c("lower.CL","upper.CL")])

pairwise.uamps.hem <- pairwise.uamps.hem[c("contrast","estimate","lower.CL",
                                           "upper.CL","t.ratio","p.value")]

pairwise.uamps.hem[,"d"] <- pairwise.uamps.hem$estimate/sigma(uamp7)

pairwise.uamps.hem[,2:7] <- round(pairwise.uamps.hem[,2:7],2)

colnames(pairwise.uamps.hem) <- c("contrast","estimate",
                                  "CI 2.5%","CI 97.5%","t","p","d")
kable(pairwise.uamps.hem)

# make table  
write.table(pairwise.uamps.hem,file="/Volumes/Projects/2017-0121-CCMusic/analysis_MM/feature_hem_interaction_amp.csv",
            sep=",",row.names = FALSE,quote=FALSE) # Export to a table


# Table 5a: pairwise contrasts of mean amplitudes for hemisphere for features / right-left interaction

confint.uamps.hem <- as.data.frame(confint(lsmeans(uamp7,pairwise~hem|dev,
                                                   adjust = "Bonferroni")$contrast))

pairwise.uamps.hem <-cbind(as.data.frame(lsmeans(uamp7, pairwise~hem|dev,
                                                 adjust="bonferroni")$contrasts),
                           confint.uamps.hem[c("lower.CL","upper.CL")])

pairwise.uamps.hem <- pairwise.uamps.hem[c("contrast","dev","estimate","lower.CL",
                                           "upper.CL","t.ratio","p.value")]

pairwise.uamps.hem[,"d"] <- pairwise.uamps.hem$estimate/sigma(uamp7)

pairwise.uamps.hem[,3:8] <- round(pairwise.uamps.hem[,3:8],2)

colnames(pairwise.uamps.hem) <- c("contrast","feature","estimate",
                                  "CI 2.5%","CI 97.5%","t","p","d")
kable(pairwise.uamps.hem)

# TABLE  
write.table(pairwise.uamps.hem,file="/Volumes/Projects/2017-0121-CCMusic/analysis_MM/feature_hem_interaction_amp.csv",
            sep=",",row.names = FALSE,quote=FALSE) # Export to a table


#########################################################
#   AMPLITUDE FIGURES for uncertainty
#########################################################

uamps <- ggplot(d_long,aes(x=cond, y=amp)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(alpha = 0.6,color = 'black',size = 0.4) +
  geom_line(aes(group = ID), alpha = 0.5, size = 0.1) +
  geom_boxplot(aes(fill = cond),alpha = 0.7,fatten = 0.9, lwd = 0.1,
               color = 'black', width = 0.15, outlier.size = 0.4) +
  geom_violin(aes(color = cond),color = 'black',
              alpha = 0.2,trim = FALSE, size = 0.15) + 
  scale_fill_manual(values = c('red','blue')) +
  facet_grid(dev~hem) +
  xlab('uncertainty') +
  ylab('mean amplitude (fT)') +
  theme_bw() +
  theme(legend.position = "none"); 

uamps

#ylab('mean AMPLITUDE (\U1D707V)') +

#############################################################################################
## LATENCIES
#############################################################################################

# 1. fit mixed effect models

ulat0 <- lmer(lat~1+(1|ID), data = d_long, REML = FALSE) # control=lmerControl(optimizer="bobyqa")) # model without any effects; null-model
ulat1 <- lmer(lat~cond+(1|ID), data = d_long, REML = FALSE) # model with effect of condition; + random effect; random intercept change the basline per subject to mimic variance in a population
ulat2 <- lmer(lat~cond+dev + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
ulat3 <- lmer(lat~cond+dev+hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
ulat4 <- lmer(lat~cond+dev+hem + cond:dev + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
ulat5 <- lmer(lat~cond+dev+hem + cond:dev + cond:hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
ulat6 <- lmer(lat~cond+dev+hem + cond*dev + cond*hem + dev*hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond
ulat7 <- lmer(lat~cond+dev+hem + cond:dev + cond:hem + dev:hem + cond:dev:hem + (1|ID), data = d_long, REML = FALSE) # random effects in brackets; across subjects (ID!); random intercept = 1; random slope = cond

# Likelihood ratio test and store models:

ulats.test <- anova(ulat0,ulat1,ulat2,ulat3,ulat4,ulat5,ulat6,ulat7); 
ulats.test

ulats.models <- list(ulat0,ulat1,ulat2,ulat3,ulat4,ulat5,ulat6,ulat7);

#################################################################################

# PAIRWISE CONTRASTS

# TABLE 3: pairwise comparisons with three factors

confint.ulats3 <- as.data.frame(confint(lsmeans(ulat7,pairwise~cond|dev|hem,
                                                adjust = "Bonferroni")$contrast))

pairwise.ulats3 <-cbind(as.data.frame(lsmeans(ulat7, pairwise~cond|dev|hem,
                                              adjust="bonferroni")$contrasts),
                        confint.ulats3[c("lower.CL","upper.CL")])

pairwise.ulats3 <- pairwise.ulats3[c("contrast","dev","hem","estimate","lower.CL",
                                     "upper.CL","t.ratio","p.value")]

pairwise.ulats3[,"Cohen's d"] <- pairwise.ulats3$estimate/sigma(ulat7)
pairwise.ulats3[,4:8] <- round(pairwise.ulats3[,4:8],2)
colnames(pairwise.ulats3) <- c("contrast","deviant","hemisphere","estimate",
                               "CI 2.5%","CI 97.5%","t","p","Cohen's d")
kable(pairwise.ulats3)

# make table
write.table(pairwise.ulats3,file= "/Volumes/Projects/2017-0121-CCMusic/analysis_MM/pairwise_cond_ulats7.csv",
            sep=",",row.names = FALSE,quote=FALSE)

#################################################################################

# TABLE 4b:  

confint.ulats2 <- as.data.frame(confint(lsmeans(ulat7,pairwise~dev,
                                                adjust = "Bonferroni")$contrast))

pairwise.ulats2 <-cbind(as.data.frame(lsmeans(ulat7, pairwise~dev,
                                              adjust="bonferroni")$contrasts),
                        confint.ulats2[c("lower.CL","upper.CL")])

pairwise.ulats2 <- pairwise.ulats2[c("contrast","estimate","lower.CL",
                                     "upper.CL","t.ratio","p.value")]

pairwise.ulats2[,"d"] <- pairwise.ulats2$estimate/sigma(ulat7)
pairwise.ulats2[,2:7] <- round(pairwise.ulats2[,2:7],2)
colnames(pairwise.ulats2) <- c("contrast","deviant","estimate",
                               "CI 2.5%","CI 97.5%","t","p","Cohen's d")
kable(pairwise.ulats2)

##########

# TABLE 5b

confint.ulats.hem <- as.data.frame(confint(lsmeans(ulat7,pairwise~hem|dev,
                                                   adjust = "Bonferroni")$contrast))

pairwise.ulats.hem <-cbind(as.data.frame(lsmeans(ulat7, pairwise~hem|dev,
                                                 adjust="bonferroni")$contrasts),
                           confint.ulats.hem[c("lower.CL","upper.CL")])

pairwise.ulats.hem <- pairwise.ulats.hem[c("contrast","dev","estimate","lower.CL",
                                           "upper.CL","t.ratio","p.value")]

pairwise.ulats.hem[,"d"] <- pairwise.ulats.hem$estimate/sigma(uamp7)

pairwise.ulats.hem[,3:8] <- round(pairwise.ulats.hem[,3:8],2)

colnames(pairwise.ulats.hem) <- c("contrast","feature","estimate",
                                  "CI 2.5%","CI 97.5%","t","p","d")
kable(pairwise.ulats.hem)


############


# FIGURE for uncertainty - latency analyses

ulats <- ggplot(d_long,aes(x=cond, y=lat)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(alpha = 0.6,color = 'black',size = 0.4) +
  geom_line(aes(group = ID), alpha = 0.5, size = 0.1) +
  geom_boxplot(aes(fill = cond),alpha = 0.7,fatten = 0.9, lwd = 0.1,
               color = 'black', width = 0.15, outlier.size = 0.4) +
  geom_violin(aes(color = cond),color = 'black',
              alpha = 0.2,trim = FALSE, size = 0.15) +
  scale_fill_manual(values = c('red','blue')) +
  facet_grid(dev~hem) +
  xlab('uncertainty') +
  ylab('peak latency (ms)') +
  theme_bw() +
  theme(legend.position = "none"); 

ulats

###############################################
# Make joint uncertainty reports:
###############################################

# mixed-effect models

# amplitudes
#uncertainty.report <- data.frame('model' = rownames(uamps.test))
#uncertainty.report[2:nrow(uncertainty.report),'null'] <- uncertainty.report[1:nrow(uncertainty.report)-1,'model']

#uncertainty.report <- cbind(uncertainty.report,round(uamps.test[,c('AIC','Chisq','Pr(>Chisq)')],2),round(ulats.test[,c('AIC','Chisq','Pr(>Chisq)')],2))
#colnames(uncertainty.report) <- c('model','null','AIC','X2','p','AIC','X2','p')

#write.table(uncertainty.report,file="/Volumes/Projects/2017-0121-CCMusic/analysis_MM/uncertainty_report.csv",
#            sep= ",",row.names = FALSE,quote=FALSE)

# pairwise contrasts:

#uncertainty.pw <- cbind(pairwise.uamps.hem,pairwise.ulats.hem[,c(2:ncol(pairwise.ulats.hem))])

#write.table(uncertainty.pw,file="/Volumes/Projects/2017-0121-CCMusic/analysis_MM/feat_hem_interaction_all.csv",
#            sep=",",row.names = FALSE,quote=FALSE)

#write.table(uncertainty.pw,file="/Users/iris.mencke/Documents/feat_hem_interaction_all.csv",
#            sep=",",row.names = FALSE,quote=FALSE)



# JOINT FIGURE

uplots <- arrangeGrob(uamps,ulats,ncol=2);
plot(uplots)

ggsave("/Volumes/Projects/2017-0121-CCMusic/analysis_MM/uncertainty.pdf", plot=uplots,width = 180, height = 190, units = 'mm', dpi = 300)
ggsave("/Volumes/Projects/2017-0121-CCMusic/analysis_MM/uncertainty.png", plot=uplots,width = 180, height = 190, units = 'mm', dpi = 300)

