rm(list=ls())
N <- 10# polpulation size
sl <- 4 # sample size

rds.sample <- data.frame(MyUniID = seq(1, by = 1, len = sl) ,
                         NS1 =  seq(2, by = 0, len = sl),
                         refCoupNum  = seq(0, by = 1, len = sl), 
                         coup1 = seq(1, by = 1, len = sl),
                         coup2 = seq(sl+1, by = 1, len = sl), 
                         coup3 = seq(2*sl+1, by = 1, len = sl),
                         interviewDt = seq(0, by = 0, len = sl), 
                         stringsAsFactors = TRUE)
rds.sample$NS1[1] <- 3
for (r in 2:sl)
{
  rds.sample$interviewDt[r] <- rds.sample$interviewDt[r-1] + 1/((r-1)*(N-r+1))
}


rds.object<- initializeRdsObject(rds.sample, seeds=1)
rds.object$estimates <- estimate.b.k(rds.object = rds.object )
signif(sum(rds.object$estimates$Nk.estimates), digits = 2)
getTheta(rds.object)
