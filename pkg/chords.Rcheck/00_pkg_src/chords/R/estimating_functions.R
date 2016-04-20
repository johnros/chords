estimate.b.k.2 <- function(k, A.k, B.k, n.k, n.k.count, k.ind, regularize=FALSE, silent){
  ## Initialize:
  result <- list(N.k=NA, converge=FALSE)
  N.k <- NA
  
  # Unregulirized estimator:
  if(!regularize){
    target <- function(N.k){
      const1 <- N.k - (B.k/A.k)
      # pre.const2 <- (N.k - n.k[which(k.ind)-1])
      pre.const2 <- (N.k - seq(0, n.k.count-1))
      
      # Deal with impossible N.k:
      if(any(pre.const2 < 0)) return(-Inf)
      
      const2 <- sum(1/pre.const2)
      const1 * const2 - n.k.count
    }
  } 
  # Regulirized case:
  else{
    # Construct target function:
    target <- function(N.k){
      const1 <- N.k-(B.k/A.k)
      # pre.const2 <- (N.k - n.k[which(k.ind)-1])
      pre.const2 <- (N.k - seq(0, n.k.count-1))
      
      # Deal with impossible N.k:
      if(any(pre.const2<0)) return(-Inf)
      
      const2.1 <- sum(1/pre.const2^3)
      const2.2 <- sum(1/pre.const2^2)
      const2.3 <- sum(1/pre.const2)
      const2 <- const2.3 - const2.1/const2.2
      const1*const2 - n.k.count
    }
  }
  
  # Estimate n by finding root of target:
  roots <- NULL
  .interval <- n.k.count * c(1, 10*max(1,1/A.k))
  try(roots <- uniroot(f = target, interval =.interval, 
                       extendInt = "no", maxiter = 1e4), silent = silent)
  
  # In case of convergence:
  if(length(roots)>0) {
    result$N.k <- roots$root
    result$converge <- 0 # 0 marks succeful convergence!
  } 
  # In case of non convergence:
  else {
    result$N.k <- Inf #n.k.count
    result$converge <- sign(target(.interval[2]))
  }
  
  result$N.k <- ceiling(result$N.k)
  return(result)
}
## Testing:
# save(n.k, k.ind,A.k, B.k, file='temp/testing_setup.RData')
# load(file='temp/testing_setup.RData')
# estimate.b.k.2(k = 3, A.k = A.k, B.k = B.k, n.k =n.k, k.ind = k.ind )






# estimate beta_k from sampled degrees and snowball matrix:
estimate.b.k<- function (rds.object, 
                         const=1,
                         impute.Nks='jeff',
                         silent=TRUE) {
  ### Sketch:
  # Generate estimable parameters vector.
  # Optimized parameter-wise.
  
  
  ### Verifications:
  if(length(rds.object$estimates)>0) {
    message('Overwriting existing estimates in rds.object.')  
  }
  
  ### Initialize:
  arrival.times <- formatArrivalTimes(rds.object$rds.sample$interviewDt)
  arrival.intervals <- diff(arrival.times)
  
  arrival.degree<- rds.object$rds.sample$NS1
  max.observed.degree<- max(arrival.degree)
  ## TODO: should the degree of the first (seed) be removed?
  degree.counts<- table(arrival.degree)
  
  # Sequences per degree
  I.t <- rds.object$I.t
  degree.in <- rds.object$degree.in
  degree.out <- rds.object$degree.out
  
  likelihood <- NA
  
  ### Estimate:
  Nk.estimates<- rep(9999L, max.observed.degree) 
  names(Nk.estimates)<- seq_len(max.observed.degree)
  log.bk.estiamtes<- rep(NA, max.observed.degree) 
  A.ks<- rep(NA, max.observed.degree) 
  B.ks<- rep(NA, max.observed.degree) 
  n.k.counts<- rep(NA, max.observed.degree)
  convergence<- rep(NA, max.observed.degree)
  
  A.k <- sum( head(I.t,-1) * arrival.intervals, na.rm=TRUE)
  
  uniques<- as.integer(names(degree.counts))
  Nk.estimates[-uniques]<- 0
  for(k in uniques){
    # k <- uniques[[1]]
    k.ind <- arrival.degree==k
    k.ind[1] <- FALSE # dealing with sample kickoff
    
    n.k <- cumsum((arrival.degree==k))
    #     n.k <- cumsum((degree.in==k) - (degree.out==k))
    n.k.count <- degree.counts[paste(k)]
    
    #     head(cbind(arrival.degree, n.k, I.t, arrival.intervals),20)
    #     head(cbind(arrival.degree[-1], n.k[-1], I.t[-1], arrival.intervals))
    #     head(cbind(arrival.degree[-1], head(n.k,-1), head(I.t,-1), arrival.intervals))
    B.k <- sum( head(I.t,-1) * arrival.intervals * head(n.k,-1), na.rm=TRUE)    
    
    
    .temp <- estimate.b.k.2(k=k, 
                            A.k=A.k*const, 
                            B.k=B.k*const, 
                            n.k=n.k, 
                            n.k.count= n.k.count, 
                            k.ind=k.ind, 
                            regularize=FALSE, 
                            silent = silent)    
    
    convergence[k] <- .temp$converge
    A.ks[k] <- A.k
    B.ks[k] <- B.k
    n.k.counts[k] <- n.k.count
    
    # Impute using jeffrey's prior:
    if(impute.Nks=='jeff' && .temp$converge==0){
      .temp <- estimate.b.k.2(k=k, A.k=A.k*const, B.k=B.k*const, n.k=n.k, 
                              n.k.count= n.k.count, k.ind=k.ind, 
                              regularize=TRUE)    
    }
    
    Nk.estimates[k]<-.temp$N.k
    log.bk.estiamtes[k] <- log(n.k.count) - log(.temp$N.k *  A.k - B.k)
  } # End looping over estimable degrees.  
  
  
  
  
  # Impute using my heuristic:
  if(impute.Nks=='john') {
    Nk.estimates <- imputeEstimates(Nk.estimates, n.k.counts, convergence)
  }

  
  
  likelihood.val <- likelihood(log.bk = log.bk.estiamtes, 
                           Nk.estimates = Nk.estimates, 
                           I.t = I.t, 
                           n.k.counts = n.k.counts, 
                           degree.in = degree.in , 
                           degree.out = degree.out , 
                           arrival.intervals = arrival.intervals, 
                           arrival.degree = arrival.degree)
  
  result<- list(
    call=sys.call(),
    Nk.estimates=Nk.estimates, 
    log.bk.estimates=log.bk.estiamtes,
    A.ks=A.ks,
    B.ks=B.ks,
    n.k.counts=n.k.counts,
    arrival.intervals=arrival.intervals,
    arrival.degree=arrival.degree, 
    convergence=convergence,
    likelihood=likelihood.val)
  
  return(result)    					
}






