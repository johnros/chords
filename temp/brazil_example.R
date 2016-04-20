#-------------- Preparing  Data -----------------#
rm(list=ls())
library(dplyr)
library(lubridate)
library(magrittr)

# source('temp/Uganda_example.R')


brazil.csv <- read.csv(file='~/Dropbox/RDS/Round2/Brazil/raw/multiplier.csv', stringsAsFactors=FALSE)
brazil.csv2 <- tbl_df(brazil.csv)

table(table(brazil.csv2$interview.date))



brazil.csv3 <- select(brazil.csv2, 
                      uid, netsize.1.bss, InCoupon, OutCoupon1, OutCoupon2, OutCoupon3, interview.date)
names(brazil.csv3) <- c("MyUniID", "NS1", "refCoupNum", "coup1", "coup2", "coup3", 
                        "interviewDt")

brazil.csv3$interviewDt2 <- brazil.csv3$interviewDt
brazil.csv3$interviewDt <- as.character(brazil.csv3$interviewDt) %>%  dmy %>% unclass

# write.table(brazil.csv3, file='pkg/chords/data/brazil.tab', row.names=FALSE)

#----------- Analyzing --------------#
rm(list=ls())
data(brazil)

table(table(brazil$interviewDt))

rds.object2<- initializeRdsObject(brazil)
rds.object2$estimates <- estimate.b.k(rds.object = rds.object2 )

# View estimates:
plot(rds.object2$estimates$Nk.estimates, type='h')
# Population size estimate:
sum(rds.object2$estimates$Nk.estimates)
plot(rds.object2$estimates$log.bk.estimates, type='h')

# Estiamte discoverability coefficient:
getTheta(rds.object2)$theta
thetaSmoothingNks(rds.object2)

# How many degrees were imputed?:
table(rds.object2$estimates$convergence)




