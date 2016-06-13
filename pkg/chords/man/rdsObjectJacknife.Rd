\name{rdsObjectJacknife}
\alias{rdsObjectJacknife}

\title{
  Update the sample for Jacknife estimation
}
\description{
  A utility function for Jacknifing the RDS likelihood estimator. 
  Removes an observation from the sample. 
  This has the effect of changing the sampling timestamps of all following observations.
}
\usage{
rdsObjectJacknife(rds.object, jack.ind, seeds = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds.object}{
  An RDS class object, as outputted by \code{\link{initializeRdsObject}}.
}
  \item{jack.ind}{
  The index of the observation to be removed.
}
  \item{seeds}{
  The number of initial recruiting agents.
}
}

\value{
An RDS class object, after removing an observation, updating timestamps, and reodring sample.
}

\references{
Musa, John D., and A. Iannino. "Estimating the Total Number of Software Failures Using an Exponential Model."
SIGSOFT Softw. Eng. Notes 16, no. 3 (July 1991): 80-84. doi:10.1145/127099.127123.
}


\seealso{
  \code{\link{initializeRdsObject}}
}

\examples{
\dontrun{
## Preliminaries:

library(chords)
library(magrittr)

dk <- c(2, 10) # unique degree classes
true.dks <- rep(0,max(dk)); true.dks[dk] <- dk
true.Nks <- rep(0,max(dk)); true.Nks[dk] <- 1e3
beta <- 1 
theta <-  0.1
true.log.bks <- rep(-Inf, max(dk))
true.log.bks[dk] <- theta*log(beta*dk)
sample.length <- 4e2
nsims <- 1e2

# Make samples
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
  simlist[[i]]$estimates <- estimate.b.k(rds.object = simlist[[i]], impute.Nks='john', silent=TRUE)
  llvec[i] <- simlist[[i]]$estimates$likelihood
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
  nvec[i] <- sum(simlist[[i]]$estimates$Nk.estimates)
  converged[i] <- sum(simlist[[i]]$estimates$convergence, na.rm=TRUE)
}
summary(beta0vec)
summary(nvec)
hist(nvec)
abline(v=sum(true.Nks), col='red')
abline(v=median(nvec, na.rm = TRUE), lty=2)
table(converged) # 0 for convergence. 


## Now finally to the Jacknife
nvec
rds.object <- simlist[[1]]
Nks <- rds.object$estimates$Nk.estimates
convergence <- rds.object$estimates$convergence
jack <- which(convergence!=0)
jack.Nks <- names(Nks)[jack] %>% as.numeric

jack.ind <- rds.object$degree.in %in% jack.Nks
# jack.ind <- rep(TRUE, nrow(rds.object$rds.sample))
N <- sum(jack.ind)
n.jack <- rep(NA, N-1)
for(j in seq_along(jack.ind)){
  # j <- 1
  i <- which(jack.ind)[j]
  rds.jack <- rdsObjectJacknife(rds.object, i)
  rds.jack.estim <- estimate.b.k(rds.jack)  
  n.jack[j] <- sum(rds.jack.estim$Nk.estimates)
}
n.jack

}

}