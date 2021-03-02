working.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working.dir)

d <- read.table('clean_data/dataset.csv', header = T, sep = ',')
d <- d[d$prof < 2,]
d <- d[d$yomt < 4,]
d <- d[d$somt > 9|d$somt ==0,]

dagg <- aggregate(d, by = list(d$sub), mean)

r <- round(data.frame(
  "age" = c(mean(dagg$age),sd(dagg$age)),
  "years of musical training" = c(mean(dagg$yomt),sd(dagg$yomt)),
  "age of start of musical training" = c(mean(dagg$somt),sd(dagg$somt))
),2)

table(d$sex)/48