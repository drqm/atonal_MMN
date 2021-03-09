setwd("C:/Users/au571303/Documents/projects/atonal_MMN/IDyOM/")
library(ggplot2)
#library(gridExtra)
library(patchwork)
files <- list.files('IDyOM_models_pclass')

for (fidx in 1:length(files)){
  f = files[fidx]
  file_info = strsplit(f,'-')
  cdata <- read.table(paste0('IDyOM_models_pclass/',f),header = TRUE)
  cdata <- cdata[,c('melody.id','cpitch','information.content','entropy')]
  cdata$condition <- file_info[[1]][1]
  cdata$model <- file_info[[1]][9]
  cdata$training_set <- file_info[[1]][4]
  if (fidx == 1){
    d <- cdata
  }else{
    d <- rbind(d,cdata)
  }
}

d$corpus <- ''
d$corpus[which(d$training_set == 'nil')] <- 'no corpus'
d$corpus[which(d$training_set == '30_31_32')] <- 'tonal corpus'
d$corpus[which(d$training_set == '34_35')] <- 'atonal corpus'
d$stim_type <- ifelse(d$condition == 28, 'tonal','atonal')

d_agg <- aggregate.data.frame(d, by = list(d$melody.id,d$stim_type,d$corpus),mean)
d_agg <- d_agg[c(2:7)]
colnames(d_agg)[1:2] <- c('condition','corpus')

d_agg$corpus <- factor(d_agg$corpus, levels = c("no corpus", "atonal corpus", "tonal corpus"))
eplot <- ggplot(d_agg, aes(condition,entropy)) +
  geom_jitter(width = 0.01, size = 2, alpha = 0.4, color = 'black') +
  geom_boxplot(alpha = 0.7, width = 0.4, aes(fill = condition)) +
  scale_fill_manual(values = c('blue','red')) + 
  facet_wrap(~corpus) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab('mean entropy per melody (bits)')

ic.plot <- ggplot(d_agg, aes(condition,information.content)) +
  geom_jitter(width = 0.02, size = 2, alpha = 0.5, color = 'black') +
  geom_boxplot(alpha = 0.7, width = 0.4, aes(fill = condition)) +
  scale_fill_manual(values = c('blue','red')) + 
  facet_wrap(~corpus) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab('mean IC per melody (bits)')

notes = c('C','C#','D','D#','E','F','F#','G','G#','A','A#','B')
pdists <- ggplot(d) +  
  geom_histogram(aes(x = cpitch, y = ..density.. , fill = stim_type),
                 binwidth = 1,alpha = 0.7, color ='black') + 
  scale_fill_manual(values = c('blue','red')) + 
  facet_wrap(~stim_type) + 
  xlab('pitch class') +
  ylab('probability of occurrence') +
  scale_x_continuous(breaks = 60:71, limits = c(59,72),labels = notes) +
  theme_bw() +  theme(legend.position = "none"); pdists

## Put figures together:
# labels = c('A','B'), ncol =1,nrow =2)
# 
# ic.plots <- plot_grid(eplot,ic.plot,plots[[3]],plots[[4]],
#                          labels = c('A','','B'), ncol =2,nrow =2)
# ic.plots <- arrangeGrob(arrangeGrob(eplot,ic.plot,ncol = 2),pdists,nrow = 2);plot(ic.plots)

ic.plots <- (eplot + ic.plot) / pdists + plot_annotation(tag_levels = 'a');ic.plots
ggsave("Fig2-IC_plots.png", plot=ic.plots,width = 180, height = 180, units = 'mm', dpi = 600)
ggsave("Fig2-IC_plots.pdf", plot=ic.plots,width = 180, height = 180, units = 'mm', dpi = 600)


