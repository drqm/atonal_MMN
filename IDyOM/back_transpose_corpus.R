library(ggplot2)
setwd("C:/Users/au571303/Dropbox/Doctorado/write_midi/")
source("write_midi.R")

main_dir <- "C:/Users/au571303/Documents/projects/atonal_mumufe_beh/"
IDyOM_files <- list.files(main_dir, pattern = ".dat")

os <- 'windows'
for (s in IDyOM_files){
  stim <- read.table(paste0(main_dir,s),header = TRUE)
  stim <- stim[,c('melody.id','keysig','cpitch')]
  stim[stim$keysig %in% c(0,2,4),'cpitch'] <- stim[stim$keysig %in% c(0,2,4),'cpitch'] + 1
  stim$keysig[is.na(stim$keysig)] <- 0
  stim$new_pitch <- (stim$cpitch - stim$keysig + 6) %% 12 + 60 # transposing to c and taking modulo
  plot <- ggplot(stim) +
    geom_histogram(aes(x = new_pitch, y = ..density.. ),binwidth = 1,color ='black',alpha = 0.5) #+
  #  facet_wrap(~keysig)
  print(plot)
  for (m in unique(stim$melody.id)){
    cmel <- stim[stim$melody.id==m,]
    out_dir <- paste0(main_dir,'midi_corpus_pclass/',substr(s,1,2),'/')
    dir.create(out_dir, showWarnings = FALSE)
    out_name <- paste0(out_dir,as.character(m))
    print(out_name)
    write.midi(out_name,cmel$new_pitch,60,0,os)
  }
}

# plot <- ggplot(stim) +
#   geom_histogram(aes(x = new_pitch, y = ..density.. ),binwidth = 1,color ='black',alpha = 0.5)