rm(list=ls())
#------ Uganda Data (?) -------------#

## Import RDS file with time stamps
col.classes <- c('integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'character')
f.loc <- paste( '~/Dropbox/RDS/RDS_shared/Uganda_example/rds_data_Uganda02.csv',sep='')
rds.sample<- read.csv(file = f.loc, colClasses=col.classes)
names(rds.sample)

# Mark tied arrival times:
ties <- (table(rds.sample$interviewDt)>1) %>% 
  which %>%
  names
tie.inds <- which(rds.sample$interviewDt==ties)

## Create big rds object:
rds.object<- initializeRdsObject(rds.sample, seeds = 1)
plot(rds.object$I.t, type='h', col='grey')


## Do N.k estimation:
rds.object$estimates <- estimate.b.k(rds.object = rds.object, impute.Nks = 'jeff')
sum(rds.object$estimates$Nk.estimates)
plot(table(rds.object$rds.sample$NS1))
plot(rds.object$estimates$Nk.estimates, type='h')


## Use b_k=beta*k^theta to smooth population size:
getTheta(rds.object, robust=TRUE)$theta
sum(thetaSmoothingNks(rds.object, robust=TRUE))



#----- Proceed to maximum likelihood estimation---------#
#  chords:::estimate.b.theta(rds.object)



#---- Aggregate degrees to lower variance -----#

bin <- 7
rds.object2<- initializeRdsObject(rds.sample, bin = bin)
rds.object$estimates <- estimate.b.k(rds.object = rds.object2 )
sum(rds.object$estimates$Nk.estimates)

plot(table(rds.object$rds.sample$NS1))
plot(rds.object$estimates$Nk.estimates, type='h')
getTheta(rds.object, bin=bin)




