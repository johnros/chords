rm(list=ls())


## Repeat and compute average MSE:

true.Nks <- rep(0,100); true.Nks[c(2,100)] <- 1000
theta <- 1e-1
true.log.bks <- rep(-Inf, 100)
true.beta <- 1
true.log.bks[c(2,100)] <- theta*log(c(2,100))+log(true.beta)
sample.length <- 1000L

rds.simulated.object <- makeRdsSample(
  N.k =true.Nks, 
  b.k = exp(true.log.bks),
  sample.length = sample.length)
rds.simulated.object$estimates <- estimate.b.k(rds.object = rds.simulated.object )
chords:::compareNkEstimate(rds.simulated.object$estimates$Nk.estimates, true.Nks)
sum(rds.simulated.object$estimates$Nk.estimates)
getTheta(rds.simulated.object)$theta


## Try maximum likelihood 
chords:::estimate.b.theta(rds.simulated.object)


true.Nks <- rep(0,100); true.Nks[c(2,100)] <- 1000
theta <- 1e-1
beta_0 <- 1e-3
true.log.bks <- rep(-Inf, 100)
true.log.bks[c(2,100)] <- log(beta_0) + theta*log(c(2,100))
sample.length <- 1000L
replicationss <- 100L
MSEs <- replicate(replicationss,{
  rds.simulated.object <- makeRdsSample(
    N.k =true.Nks , 
    b.k = exp(true.log.bks),
    sample.length = sample.length)
  rds.simulated.object$estimates <- estimate.b.k(rds.object = rds.simulated.object )
  chords:::compareNkEstimate(rds.simulated.object$estimates$Nk.estimates, true.Nks)
})
MSEs
plot(t(MSEs), xlim=c(0,10), ylim=c(0,10), main=paste('theta=',theta,sep=''));abline(h=1, v=1)





replicationss <- 100L
thetas <- replicate(replicationss,{
  rds.simulated.object <- makeRdsSample(
    N.k =true.Nks , 
    b.k = exp(true.log.bks),
    sample.length = sample.length)
  rds.simulated.object$estimates <- estimate.b.k(rds.object = rds.simulated.object )
  getTheta(rds.simulated.object)
})
plot(unlist(thetas[2,]), ylim=c(0,1))
abline(h=0.1); abline(h=median(unlist(thetas[2,])),lty=2)



source('temp/Uganda_example.R')
rds.object$estimates <- estimate.b.k(rds.object = rds.object, impute.Nks = FALSE )



