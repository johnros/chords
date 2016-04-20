# 
# Author: johnros
###############################################################################




generate.rds.control<- function(				
){
	return(list(NULL))
}
## Testing:
# generate.rds.control()





# Solve for a fixed Nj given theta and beta:
NjSolve <- function(sampled.degree.vector, S, Sij, j, Nj.table, beta, theta, maximal.Nj=1e6, ...){
	## Initialize:
	result<- NA
	
	const1 <- beta * j^theta
	
	target <- function(Nj){
		
		Uij <- Nj-Sij[paste(j),]
		sum( (sampled.degree.vector==j) / Uij - (sampled.degree.vector!=j) * const1 * S / (1 - const1 * S* Uij) )
	}
	
	
	Nj.bound<- min((1+beta*j^theta*S*Sij[paste(j),])/(beta*j^theta*S))
	Nj.low <- Nj.table[paste(j)]+1
	Nj.high <- floor(Nj.bound)
	
	if(Nj.high < Nj.low) return(result)
	
	try(result <- uniroot(target, interval = c(Nj.low, Nj.high),...), silent=TRUE)
	
	if(length(result)>1) return(result$root)
	
	else if ( !is.infinite(target(Nj.low)) && !is.infinite(target(Nj.high)) && target(Nj.high)*target(Nj.low)>0  ){
		if(target(Nj.low) < 0) return(Nj.low)
		if(target(Nj.low) > 0) return(Nj.high)
	}
	
}
## Testing:
# Njs<- c(100,200,200,200)
# names(Njs)<- c("10","50","100","1000")
# Njs<- as.table(Njs)
# theta<- 1.1
# beta<- 3e-9
# degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4)
# (.Nj.table <- table(degree.sampled.vec))
# matplot(t(.Sij <- make.Sij(degree.sampled.vec)), type="s")
# plot(.S <- compute.S(.Sij), type="s")
# undebug(chords:::NjSolve)
# chords:::NjSolve(degree.sampled.vec, .S, .Sij, Nj.table = .Nj.table, j=100, beta = beta, theta = theta, maximal.Nj = 1e15)
# chords:::NjSolve(degree.sampled.vec, .S, .Sij, Nj.table = .Nj.table, j=50, beta = beta, theta = theta, maximal.Nj = 1e15)






# Get a list of beta, theta, returns a *table* with the MLE of every *observed* degree:
estimateNjs <- function(x, sampled.degree.vector, S, Sij, Nj.observed.table){
	# Note: Nj.table is expected to have zero counts *removed*
	
	beta <- unlist(x['beta'])
	theta <- unlist(x['theta'])
	
	## TODO: Check if theta and beta values  
	
	
	Njs <- rep(NA, length(Nj.observed.table))
	names(Njs) <- names(Nj.observed.table)
	
	for (j in as.numeric(names(Nj.observed.table))){
    # j<- 10
		Njs[paste(j)] <- NjSolve(sampled.degree.vector, S=S, Sij=Sij, Nj.table=Nj.observed.table, j=j, beta = beta, theta = theta, maximal.Nj = 1e15)
	}
	
	
	return(list(beta=beta, theta=theta, Njs=Njs))	
}
## Testing:
# Njs<- c(100,200,200,200)
# names(Njs)<- c("10","50","100","1000")
# Njs<- as.table(Njs)
# theta<- 1.1
# beta<- 3e-9
# degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4)
# (.Nj.table <- table(degree.sampled.vec))
# matplot(t(.Sij <- make.Sij(degree.sampled.vec)), type="s")
# plot(.S <- compute.S(.Sij), type="s")
# 
# debug(chords:::estimateNjs)
# chords:::estimateNjs(list(beta=beta, theta=theta), degree.sampled.vec, .S, .Sij, .Nj.table[-1])








getBetaBound <- function(Nj.table, Sij, S, theta) {
	beta.bound.j <- function(j){
		Uij <- Nj.table[paste(j)] - Sij[paste(j),]
		beta.high<-   min( -1 * (log(S) + theta * log(j) + log(Uij) ), na.rm=TRUE)
		return(beta.high)
	}
	beta.bound <- min(sapply(as.numeric(names(Nj.table)), beta.bound.j), na.rm=TRUE)
	return(beta.bound)
}
## Testing:
#Njs<- c(100,200,200,200)
#names(Njs)<- c("10","50","100","1000")
#Njs<- as.table(Njs)
#theta<- 1.1
#beta<- 3e-9
#degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4)
#(.Nj.table <- table(degree.sampled.vec))
#matplot(t(.Sij <- make.Sij(degree.sampled.vec)), type="s")
#plot(.S <- compute.S(.Sij), type="s")

# chords:::getBetaBound(.Nj.table[-1], .Sij, .S, theta)






makeBetaThetaGrid <- function(thetas, beta.length, sampled.degree.vector, S, Sij, Nj.observed.table){
	result <- list()
	
	for(i in 1:length(thetas)){
		.theta <- thetas[i]
		beta.high.log <- getBetaBound(Nj.table = Nj.observed.table[-1], Sij, S, .theta)
		log.betas <- seq(from=beta.high.log, by=-0.01, length=beta.length)				
		
		for(k in 1:beta.length){
			result <- c(result, list(c(theta=thetas[i], beta=exp(log.betas[k]))))			
		}
	}
	return(result)	
}
## Testing:
#Njs<- c(100,200,200,200)
#names(Njs)<- c("10","50","100","1000")
#Njs<- as.table(Njs)
#theta<- 1.1
#beta<- 3e-9
#degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4)
#(.Nj.table <- table(degree.sampled.vec))
#matplot(t(.Sij <- make.Sij(degree.sampled.vec)), type="s")
#plot(.S <- compute.S(.Sij), type="s")

#chords:::makeBetaThetaGrid(thetas=seq(-1,1,length=10), beta.length=10, degree.sampled.vec, .S, .Sij, .Nj.table[-1])








## Estimate using a grid over theta and beta and find root of derivative for Nj:
estimate.rds<- function (sampled.degree.vector, Sij, initial.thetas, beta.grid.length, 
		arc=FALSE, control=generate.rds.control()) {
	
	## Initializing:
	the.call<- sys.call()
	
	
	max.observed.degree<- max(sampled.degree.vector)
	N.j<- rep(0, max.observed.degree)
	# Look for the degrees in the data with non trivial estimates:
	Nj.observed.table<-table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js<- as.numeric(names(Nj.observed.table))
	maximal.degree.count<- max(Nj.observed.table)
	final.result<- NA
	
	
	# Compute the size of the snowball along the sample: (aka I_t)
	S<- compute.S(Sij)
	
	### Generate theta and beta grid:
	# Fills and formats initialization values as required by optim.wrap()
	# Assumes initial.values is a list of lists. Each containing a theta, beta and an Nj vector.	
		 
	beta.theta.grid <- makeBetaThetaGrid(
			thetas=initial.thetas, 
			beta.length=beta.grid.length,
			sampled.degree.vector = sampled.degree.vector,
			S=S, 
			Sij= Sij, 
			Nj.observed.table=Nj.observed.table	)
	
	
	
	#### Estimate Nj given theta and beta
	estimateNjs.local <- function(x) try(estimateNjs(x, sampled.degree.vector = sampled.degree.vector, S=S, Sij=Sij, Nj.observed.table = Nj.observed.table), silent=TRUE)
	beta.theta.Nj.grid <- lapply(beta.theta.grid, estimateNjs.local)
	
	

	# Wrap the likelihood function. 
	# Implements constraints on parameters by taking real valued (canonical) parameters and remapping them.
	likelihood.wrap<- function(x){
		# Initialize:
		beta <- x$beta
		theta <- x$theta
		Nj.table <- x$Nj
		likelihood.result<- -1e15
		
		if(is.na(beta) || is.na(theta) || any(is.na(Nj.table))) return(likelihood.result)
		
		
		N.j<- rep(0, max.observed.degree)
		N.j[Observed.js]<- Nj.table[paste(Observed.js)] # fill non trivial Nj estimates.
		
		
		# Checking estimates are within allowed range:
		bad.Nj<- any( round(N.j[Observed.js],10) < round(Nj.observed.table,10) ) # Estimated Nj no smaller than observed
		
		log.beta.bound <- getBetaBound(Nj.observed.table, Sij, S, theta)
		bad.beta<- beta > exp(log.beta.bound)  # Defined beta values
		if(isTRUE( bad.Nj || bad.beta)) return(likelihood.result) 
		
		# Compute likelihood:
		if(is.numeric(beta) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(beta) ) {
			result<- .C("likelihood", 
					sample=as.integer(sampled.degree.vector), 
					Sij=as.integer(as.matrix(Sij)),
					S=as.integer(S),				  
					beta=as.numeric(beta), 
					theta=as.numeric(theta), 
					Nj=as.numeric(N.j), 
					observed_degrees=as.integer(rownames(Sij)),
					n=as.integer(length(sampled.degree.vector)), 
					N=as.integer(length(N.j)),
					N_observed=as.integer(nrow(Sij)),				  
					arc=arc,
					result=as.double(0))		  
			
		}	
		if(!is.infinite(result$result)) likelihood.result<-result$result
		return(likelihood.result)		
	}	
	
	
	#### Compute likelihood over all grid
	likelihoods <- lapply(beta.theta.Nj.grid, likelihood.wrap)
	best.index <- which.max(unlist(likelihoods))
	
	
	
	#### Prepare and return result:		
	final.result<- list( 
			beta=beta.theta.Nj.grid[[best.index]]$beta, 
			theta=beta.theta.Nj.grid[[best.index]]$theta,
			Nj=beta.theta.Nj.grid[[best.index]]$Nj,
			likelihood.optimum=likelihoods[[best.index]], 
			call=the.call)
	
	return(final.result)							
}

### Testing: 
#require(chords)
#Njs<- c(100,1000,200,200)
#names(Njs)<- c("10","20","100","1000")
#Njs<- as.table(Njs)
#theta<- 1.5
#beta<- 3e-9
#degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3)
#(.Nj.table <- table(degree.sampled.vec))
#matplot(t(.Sij <- make.Sij(degree.sampled.vec)), type="s")
#plot(.S <- compute.S(.Sij), type="s")
# estimate.rds(degree.sampled.vec, .Sij, initial.thetas=seq(-2,2,length=40), beta.grid.length=20)

# Testing with White data:
#source('~/Dropbox/Yakir/White/getData.R')
#data.Sjt <- make.Sij(data.degree)
#estimate <- estimate.rds(sampled.degree.vector=data.degree, 
#		Sij = data.Sjt, 
#		initial.thetas=seq(-2,2,length=20),
#		beta.grid.length=10,                            
#		control = generate.rds.control())









