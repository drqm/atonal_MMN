working.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working.dir)

files <- list.files('raw_data/')

# load all pilot files:
d <- data.frame()
scount <- 0
for (f in files){
  scount <- scount + 1
  cdata <- read.table(paste0('raw_data/',f),sep = ',', header = TRUE)
  cdata$participant. <- scount
  cdata[,c('key_resp_8.keys', 'key_resp_7.keys')] <- cdata[1,c('key_resp_8.keys', 'key_resp_7.keys')]
  cdata <- cdata[-c(1:5),]
  if (nrow(d) < 1){
    d <- cdata
  }else{
    d <- rbind(d,cdata[,colnames(d)])
  }
}

d <- d[,c('participant.','age.','sex..m.f..','years.of.musical.training.',
          'I.started.musical.training.at.the.age.of..',
          'key_resp_8.keys', 'key_resp_7.keys','trials.thisN','trials.thisIndex','mels',
          'key_resp_2.keys','key_resp_2.rt','key_resp_3.keys',
          'key_resp_3.rt')]

colnames(d) <- c('sub','age','sex','yomt','somt','prof','style','trial_No',
                 'trial_idx','mel','resp','resp_rt','conf','conf_rt')

d$mel <- gsub('mels/','',d$mel)
rownames(d) <- 1:nrow(d)

d_splitted <- t(as.data.frame(strsplit(as.character(d$mel),"_")))
d[,c('cond','mel','deviance')] <- d_splitted
d$deviance <- gsub('.wav','',d$deviance)
d$dev <- ifelse(grepl( "dev",d$deviance),1,2)
d$acc <-  as.numeric(d$resp == d$dev)

yomt_expressions <- c("0 not even one day",'keine','?')
somt_expressions <- c("-",'/','keine','N/A','n','never','no musical training','?')

d$yomt[d$yomt %in% yomt_expressions] <- 0 
d$somt[d$somt %in% somt_expressions] <- 0 
d$somt[d$somt == '7-9'] <- 8
d$sex[d$sex == 'female'] <- 'f'
d$sex[d$sex == 'M'] <- 'm'
d$sex <- as.character(d$sex)

write.table(d,'clean_data/dataset.csv',sep = ',',row.names = F)