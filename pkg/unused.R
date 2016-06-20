



likelihoodTheta <- function(
  beta, theta, Nk.estimates, I.t, 
  n.k.counts, degree.in, degree.out, 
  arrival.intervals, arrival.degree){

  ### Initialization:
  uniques <- which(!is.na(n.k.counts))
  n.k.t <- makeNKT(uniques, degree.in, degree.out)
  betas <- beta * uniques ^ theta
  
  
  ## Computation
  result <- 0.
  for(i in seq_along(arrival.degree)){
    if(i==1) next()
  
    for(j in seq_along(uniques)){ 
      k <- uniques[[j]]
      lambda <-  betas[j] * (Nk.estimates[k] - n.k.t[j,i-1]) * I.t[i-1]
      lambda <- max(lambda, .Machine$double.eps) 
      
      A <- ifelse(arrival.degree[i]==k, log(lambda), 0) 
      B <- lambda * arrival.intervals[i-1] 
      result <- result + A - B 
    }
  }
  
  return(result)
}
## Testing:
## TODO: adapting to Simon's impementation





wrap.likelihood <- function(beta, theta, N.k, rds.object){
  I.t <- rds.object$I.t
  n.k.counts <- rds.object$estimates$n.k.counts
  degree.in <- rds.object$degree.in
  degree.out <- rds.object$degree.out
  arrival.intervals <- rds.object$estimates$arrival.intervals
  arrival.degree <- rds.object$estimates$arrival.degree
  
  likelihoodTheta(beta, theta, N.k, I.t, 
             n.k.counts, degree.in, degree.out, arrival.intervals, arrival.degree)
  
}
##Testing:
# likelihoodTheta <- chords:::likelihoodTheta
# chords:::wrap.likelihood(beta, theta, rds.object$estimates$Nk.estimates, rds.object )


estimate.b.theta <-function(rds.object,...){
  
  ## Initialize:  
  theta_0 <- getTheta(rds.object)
  beta <- exp(theta_0$log.beta_0)
  theta <- theta_0$theta
  N.k <- rds.object$estimates$Nk.estimates
  N.k.ind <- N.k!=0
  
  beta.f <- function(beta) log(beta)
  beta.inv.f <- function(beta.converted) exp(beta.converted )
  
  theta.f <- function(theta) (theta)
  theta.inv.f <- function(theta.converted) (theta.converted)
  
  N.k.f <- function(N.k) log(N.k)
  N.k.inv.f <- function(N.k.converted) exp(N.k.converted)
  
  target <- function(x){
    result <- Inf # for a minimization problem.
    beta <- beta.inv.f(x[1]) 
    theta <- theta.inv.f(x[2])
    N.k <- rep(0, length(N.k.ind))
    N.k[N.k.ind] <- N.k.inv.f(x[-c(1,2)])
    try(result <- -wrap.likelihood(beta, theta, N.k, rds.object), silent = TRUE)
    return(result)
  }
  
  init <- c(beta.converted=beta.f(beta),
            theta.converted=theta.f(theta), 
            Nks.converted=N.k.f(N.k[N.k.ind]))
# target(init)
  optimal <- optim(par =init ,fn = target, ... )  
  
  new.beta <- beta.inv.f(optimal$par[1])
  new.theta <- theta.inv.f(optimal$par[2])
  new.N.k <- rep(0, length(N.k.ind))
  new.N.k[N.k.ind] <- N.k.inv.f(optimal$par[-c(1,2)])
  #   sum(new.N.k)
  
  result <- list(
    beta=new.beta, 
    theta=new.theta, 
    N.k=new.N.k, 
    likelihood=optimal$value,
    optim.result=optimal)
  
  return(result)
} 
## Testing:
# new.estimates <- chords:::estimate.b.theta(rds.object, control=list(maxit=1e3))
# sum(new.estimates$N.k)


## Simon's likelihood:
# dkt = degrees of each individual over time
# tt = timelist
# scale parameters
# logpar = are parameters log transformed? exclude theta
# constrain = if true, n1 and n2 represent the 'excess' over the sampled numbers
likelihood.theta.2 <- function (log.bk, Nk.estimates, I.t, 
                                n.k.counts, degree.in, degree.out, 
                                arrival.intervals, arrival.degree, const=1) {
  Sk <- Nk.estimates[Nk.estimates>0]
  dk <- unique(arrival.degree)
  ndeg <- length(dk)
  N <- sum(Nk.estimates)
  res <- 0
  nstep <- length(arrival.degree)
  # start by subtracting seed
  It <- rep(0,nstep)
  Ut <- rep(0,nstep)
  # choose first obs as seed
  seeddeg <- arrival.degree[1]
  Sk[match(seeddeg,dk)] <- Sk[match(seeddeg,dk)]-1
  It[1] <- 1
  Ut[1] <- 1
  betas <- exp(log.bk[!is.na(log.bk)])
  
  # run through timesteps and calculate loglik
  for(i in 2:nstep){
    dkrates <- betas * It[i-1] * Sk
    # total rates
    sumrates <- sum(dkrates)
    deltat <- arrival.intervals[i-1]
    res <- res + dexp(deltat, sumrates, log=TRUE)
    thisdeg <- arrival.degree[i]
    thisdegidx <- match(thisdeg,dk)
    .prob <- dkrates[thisdegidx]/sumrates
    res <- res + dbinom(1, 1, prob=.prob, log=TRUE)
    Sk[match(thisdeg,dk)] <- Sk[match(thisdeg,dk)]-1
    It[i] <- It[i-1]+1
    Ut[i] <- Ut[i-1]+1
  }
  const*res
}
## Testing:
# likelihood.theta.2(
#   log.bk = rds.simulated.object$estimates$log.bk.estimates, 
#   Nk.estimates = rds.simulated.object$estimates$Nk.estimates, 
#   I.t = rds.simulated.object$I.t, 
#   n.k.counts = rds.simulated.object$estimates$n.k.counts, 
#   degree.in = rds.simulated.object$degree.in, 
#   degree.out = rds.simulated.object$degree.out, 
#   arrival.intervals = rds.simulated.object$estimates$arrival.intervals, 
#   arrival.degree = rds.simulated.object$estimates$arrival.degree,
#   const = 1)


# generate sample:
makeWeakStrongSample <- function (means, pop.size, sample.length) {
  weights <- makeWeightMatrix(pop.size, means)
  
  ## Initiate object:
  n.ks <- rep(0, length(N.k))
  event.time <- 0
  I.t <- 1  
  rds.object.simulated <- rdsObjectConstructor(
    rds.sample = data.frame(NS1=rep(NA, sample.length), 
                            interviewDt=rep(NA, sample.length)),
    I.t = rep(NA, sample.length),
    degree.in = rep(0, sample.length),
    degree.out = rep(0, sample.length))
  
  
  ## Construct population:  
  A <- makeWeightMatrix(means, pop.size)
  
  ## Construct sample
  for(period in seq_len(sample.length)){
    ## Update sample
    ## TODO: finish sampling
    lambdas <- updateLambdaNetwork(A)
    time.to.event <- rexp(1, sum(lambda.k))
    event.time <- event.time + time.to.event
    event.type <- rmultinom(1, 1, prob = lambda.k)
    in.degree <-  which(as.logical(event.type))
    n.ks <- n.ks + event.type  
    I.t <- I.t+1
    
    
    ## Update RDSObject
    rds.object.simulated$rds.sample$interviewDt[period] <- event.time
    rds.object.simulated$degree.in[period] <- in.degree
    rds.object.simulated$I.t[period] <- I.t
    rds.object.simulated$rds.sample$NS1[period] <- in.degree
  }
  return(rds.object.simulated)
}
## Testing:
# rds.simulated.object <- makeRdsSample(N.k =N.k , sample.length = 900L)
# plot(rds.simulated.object$rds.object$rds.sample$interviewDt)

# Sketch:
## Get population A, and indices of sampled nodes
## Compute sum of weights from snowball to population.
## return vector of weights
updateLambdaNetwork <- function(A, lambdas, indices.snowball, indices.t){
  out <- A[indices.snowball,indices.t]
  enter <- 
    A[-indices.snowball, indices.t] %>% 
    result <-   lambdas + enter -out  
}
## Testing
# pop.size <- 2e3
# means <- c(weak=10, strong=2)
# A <- chords:::makeWeightMatrix(means, pop.size)
# indices.snowball <- sample(pop.size, 200)
# indices.t <- sample(pop.size, 1)
