library(magrittr)

## Simulated data example:
dk <- c(2, 1e1) # unique degree classes
true.dks <- rep(0,max(dk)); true.dks[dk] <- dk
true.Nks <- rep(0,max(dk)); true.Nks[dk] <- 1e3
beta <- 1 #5e-6
theta <-  0.1
true.log.bks <- rep(-Inf, max(dk))
true.log.bks[dk] <- theta*log(beta*dk)
sample.length <- 4e2
nsims <- 1e2

simlist <- list()
for(i in 1:nsims){
  simlist[[i]] <- makeRdsSample(
    N.k =true.Nks , 
    b.k = exp(true.log.bks),
    sample.length = sample.length)
}


# Estimate betas and theta with chords:
llvec <- rep(NA,nsims)
bklist <- list()
for(i in 1:nsims){
  # i <- 2
  simlist[[i]] <- Estimate.b.k(rds.object = simlist[[i]])
  # llvec[i] <- simlist[[i]]$estimates$likelihood
  bklist[[i]] <- simlist[[i]]$estimates$log.bk.estimates
}
b1vec <- bklist %>% lapply(.,"[",2) %>% unlist
b2vec <- bklist %>% lapply(.,"[",10) %>% unlist

hist(b1vec)
abline(v=true.log.bks[2])
hist(b2vec)
abline(v=true.log.bks[10])

beta0vec <- rep(-Inf,nsims)
thetavec <- rep(-Inf,nsims)
nvec <- rep(-Inf,nsims)
converged <- rep(9999,nsims)

for(i in 1:nsims){
  # i <- 2
  nvec[i] <- sum(simlist[[i]]$estimates$Nk.estimates)
  converged[i] <- sum(simlist[[i]]$estimates$convergence, na.rm=TRUE)
  # tfit <- getTheta(simlist[[i]])
  # beta0vec[i] <- tfit$log.beta_0
  # thetavec[i] <- tfit$theta
}
summary(beta0vec)
summary(nvec)
# summary(thetavec)
# hist(thetavec)
# abline(v=theta)
hist(nvec)
abline(v=sum(true.Nks), col='red')
abline(v=median(nvec, na.rm = TRUE), lty=2)
table(converged)

# Try various re-estimatinons:
rds.object2 <- simlist[[which(is.infinite(nvec))[1]]]

rds.object <- Estimate.b.k(rds.object = rds.object2 )
see(rds.object)
rds.object$estimates$Nk.estimates

rds.object.5 <- Estimate.b.k(rds.object = rds.object, type='rescaling')
see(rds.object.5) # will not work. less than 2 converging estimates.
rds.object.5$estimates$Nk.estimates

rds.object.6 <- Estimate.b.k(rds.object = rds.object, type='parametric')
see(rds.object.6) # will not work. less than 2 converging estimates.
rds.object.6$estimates$Nk.estimates


jack.control <- makeJackControl(3, 1e2)
rds.object.7 <- Estimate.b.k(rds.object = rds.object, 
                             type='leave-d-out', 
                             jack.control = jack.control)
see(rds.object.7)
rds.object.7$estimates$Nk.estimates


rds.object.8 <- Estimate.b.k(rds.object = rds.object, type='integrated')
see(rds.object.8)
rds.object.8$estimates$Nk.estimates


rds.object.9 <- Estimate.b.k(rds.object = rds.object, type='jeffreys')
see(rds.object.9)
rds.object.9$estimates$Nk.estimates
