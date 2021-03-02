setwd("C:/Users/au571303/Dropbox/Doctorado/write_midi/")
source("write_midi.R")
library(ggplot2)
main_dir <- "C:/Users/au571303/Documents/projects/atonal_mumufe_beh/"
os <- 'windows'
midi_constant <- 58

tonal <- list(c(2,  6,  5,  6,  9,  8,  9, 11,  9, 18, 16, 14, 13, 11,  9,  7,
   6, 14,  9,  6,  5,  6,  9,  6,  2,  1,  2,  4,  2,  6,  9,  6),
  c(21, 18, 19, 16, 14,  9, 14, 16, 18, 16, 18, 19, 21, 26, 21, 19,
   21, 18, 19, 16, 14, 18, 14,  9, 14, 16, 14, 13, 14,  6,  9,  6),
  c(2,  6,  9, 14, 13, 14, 11, 13,  9, 14, 18, 21, 19, 21, 18, 19,
    16, 14, 13, 11,  9,  7,  6,  4,  2,  6,  4,  1,  2,  6,  9,  6),
  c(2, 14, 12, 10,  9,  7,  5,  4,  5, 17, 16, 14, 12, 10,  9,  8,
    9, 21, 19, 17, 16, 14, 12, 10,  9,  7,  5,  4,  2,  5,  9,  5),
  c(2,  1,  2,  4,  5,  4,  5,  7,  9,  7,  9, 11, 13, 14, 16, 17,
    19, 17, 16, 14, 12, 10,  9,  8,  9,  7,  5,  4,  2,  5,  9,  5),
  c( 2,  5,  9,  7,  5, 14,  9, 17, 16, 14, 12, 10,  9,  7,  5,  4,
    2,  1,  2,  4,  5,  4,  5,  7,  9,  7,  9, 13, 14,  5,  9,  5))

atonal <- list(c(14,  8, 10, 12,  4,  6, 11,  7,  5, 13,  9, 15, 13,  5,  7, 11,
                6,  4, 12, 10,  8, 14,  6,  5, 12, 10,  4,  9, 11,  3,  7,  1),
               c(16, 12, 11, 15, 14, 13,  7,  5,  8,  6, 10,  9,  6,  8,  5,  7,
                13, 14, 15, 11, 12, 16,  9,  5,  6,  7, 13, 15, 12, 14, 10, 11),
               c(14, 13,  8,  7,  6,  5,  3, 12, 11,  9,  4, 10,  9, 11, 12,  3,
                5,  6,  7,  8, 13, 14,  8,  9, 10, 11, 13,  4,  5,  7, 12,  6),
               c(7,  4,  3,  9, 12,  5,  6, 11, 10,  1,  8,  2,  1, 10, 11,  6,
                 5, 12,  9,  3,  4,  7, 11, 17, 16,  9,  8, 15, 16, 13, 18, 14),
               c(15, 21, 19, 22, 20, 16, 23, 24, 18, 17, 14, 13, 17, 18, 24, 23,
                16, 20, 22, 19, 21, 15, 23, 20, 22, 26, 19, 18, 25, 26, 18, 20),
               c(9,  4, 11, 10,  8,  5, 12, 14,  6,  7, 13, 15,  8,  6, 14, 12,
                 5,  8, 10, 11,  4,  9,  7,  8, 10, 13,  6,  4, 12, 11,  5,  3))

all_mels <- list(tonal = tonal,atonal = atonal)

for (c in names(all_mels)){
  cond <- all_mels[[c]]
  out_dir <- paste0(main_dir,'midi_mels_pclass/',c,'/')
  dir.create(out_dir, showWarnings = FALSE)
  for (m in 1:length(cond)){
    cmel <- (cond[[m]] + midi_constant) %% 12 + 60
    out_name <- paste0(out_dir,as.character(m))
    print(c)
    print(m)
    print(out_name)
    write.midi(out_name,cmel,60,0,os)    
  }
}

tpitches <- data.frame('pitch'= (unlist(tonal) + midi_constant)%%12 + 60,
                       'condition' = 'tonal')
apitches <- data.frame('pitch'= (unlist(atonal) + midi_constant)%%12 + 60,
                       'condition' = 'atonal')

d <- rbind(tpitches,apitches)

ggplot(d) + 
  geom_histogram(aes(x = pitch, y = ..density.. ),binwidth = 1,color ='black',alpha = 0.5) +
  facet_wrap(~condition)

