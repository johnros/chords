# 
# Author: johnros
###############################################################################



inv.transform.theta<- function(canonical.theta, ...){
	c(canonical.theta, use.names=FALSE)
}


transform.theta<- function(theta, ...){
	c(theta, use.names=FALSE)
}




# Transforms beta back and forth to canonical (real) scale:
generate.beta.scale<- function(scale=10e10) return(scale)

transform.beta <- function(beta, scale=generate.beta.scale()){
	canonical.beta<- log(beta*scale)
	return(c(canonical.beta, use.names=FALSE))
}

inv.transform.beta<- function(canonical.beta, scale=generate.beta.scale()){
	beta<- exp(canonical.beta)/scale
	return(c(beta, use.names=FALSE))
}




# Transforms Nj back and forth to canonical (real) scale:
generate.Nj.scale <- function(scale=1) return(scale)

transform.Nj <- function(Nj, scale=generate.Nj.scale()){
	canonical.Nj<- log(Nj*scale)
	return(canonical.Nj)
}

inv.transform.Nj<- function(canonical.Nj, scale=generate.Nj.scale()){
	Nj<- exp(canonical.Nj)/scale
	return(Nj)
}
## Testing:
#inv.transform.Nj(transform.Nj(100))



#' Compares the output of the optimizaton for different initialization values.
#' @param estimation An rds2 estimation object
#' @return Used only for printing
#' @author johnros
#' @export
comparison <- function (estimation) {
	non.na<- !is.na(estimation)
	print(signif(cbind(sapply(estimation[non.na], function(x) return(x$Nj)))))
	print(sapply(estimation[non.na], function(x) sum(x$Nj)))
	print(sapply(estimation[non.na], function(x) x$theta))
	print(sapply(estimation[non.na], function(x) x[[4]]$value))
}



generate.rds.control<- function(
		maxit=500, 
		method="Nelder-Mead",
		Nj.inflations=c(1,2,10),
		beta.inflations= exp(-seq(2,0, by=-0.02)),
		fnscale=-1000,
		beta.scale=10e5
){
	return(list(maxit=maxit, method=method, Nj.inflations=Nj.inflations, beta.inflations=beta.inflations, fnscale=fnscale, beta.scale=beta.scale))
}
## Testing:
#generate.rds.control()







# Generate a list of beta values for each theta and Njs
makeBetas<- function(theta, Nj, Js, beta.inflations){
	result<-list()
	
	beta.bound<- BetaBound(Js, Nj, theta)
	for (inflation in beta.inflations){				
		result<- c(result, list(c(beta=beta.bound*inflation)))			
	}	
	return(result)
}
## Testing:
#makeBetas(theta=1,Nj=c(100,100,100), Js=c(10,50,100), beta.inflations=c(1,1e4,1e8))











makeNjs <- function(sampled.degree.vector, Nj.inflations){
	# Empirical counts
	Observed.Njs <- table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js <- as.numeric(names(Observed.Njs))
	
	# Inflate empirical counts (arbitrary)
	result<- list()
	for (inflation in Nj.inflations){
		Njs<-  ceiling(Observed.Njs*inflation)
		names(Njs)<- Observed.js
		Njs.list<- list(Njs)
		result<- c(result, Nj=Njs.list)
	}	
	# Return list of options:
	return(result)
}
### Testing:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1
#beta<- 2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#makeNjs(degree.sampled.vec, c(1,1.5))










# Generate a list of initialization values:
makeCanonicalInitialization <- function(initial.values, sampled.degree.vector, Nj.inflations, beta.inflations, scale) {
	result<- list()
	
	## Initialize with a thetas only (realistic scenario)	
	if(length(initial.values[[1]])==1L){
		stopifnot(max(sapply(initial.values, length))==min(sapply(initial.values, length)))
		Njs<- makeNjs(sampled.degree.vector, Nj.inflations)		
		# for each theta, compute a grid of betas and Njs					
		for(i in 1:length(initial.values)){
			theta<- initial.values[[i]]$theta										
			for(Nj in Njs){
				betas<- makeBetas(theta=theta, Nj=Nj, Js=as.numeric(names(Nj)), beta.inflations=beta.inflations)
				for(beta in betas){						
					# return in list-of-lists format in canonical scale!
					theta.beta.Nj<- list(
							canonical.theta=transform.theta(theta), 
							canonical.beta=transform.beta(beta,scale), 
							canonical.Nj=transform.Nj(Nj))					
					result<- c(result, list(theta.beta.Nj))						
				}					
			}			
		}			
	}
	
	
	## Initialize with theta and beta (for grid searching)		
	else if(length(initial.values[[1]])==2L){
		stopifnot(max(sapply(initial.values, length))==min(sapply(initial.values, length)))
		Njs<- makeNjs(sampled.degree.vector, Nj.inflations)	
		
		# for each theta and beta, compute a grid of betas and Njs				
		for(i in 1:length(initial.values)){
			theta<- initial.values[[i]]$theta
			beta<- initial.values[[i]]$beta
			for(Nj in Njs){					
				# return in list-of-lists format in canonical scale!
				theta.beta.Nj<- list(
						canonical.theta=transform.theta(theta), 
						canonical.beta=transform.beta(beta, scale), 
						canonical.Nj=transform.Nj(Nj))					
				result<- c(result, list(theta.beta.Nj))						
			}			
		}			
	}
	
	## Intialize with all values (for simulation verification)
	else if(length(initial.values[[1]])==3L){
		stopifnot(max(sapply(initial.values, length))==min(sapply(initial.values, length)))
		# for each theta and beta, compute a grid of betas and Njs				
		for(i in 1:length(initial.values)){
			theta<- initial.values[[i]]$theta
			beta<- initial.values[[i]]$beta
			Nj<- initial.values[[i]]$Njs				
			# return in list-of-lists format in canonical scale!
			theta.beta.Nj<- list(
					canonical.theta=transform.theta(theta), 
					canonical.beta=transform.beta(beta, scale), 
					canonical.Nj=transform.Nj(Nj))					
			result<- c(result, list(theta.beta.Nj))						
		}			
	}
	
	else stop("Bad initialization values")
	
	return(result)
}
### Testing:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1
#beta<- 2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#makeCanonicalInitialization(list(list(theta=10), list(theta=20)), degree.sampled.vec, c(1,2,3), c(1,2), 1)[[2]]
#
#makeCanonicalInitialization(list(list(theta=10)), degree.sampled.vec, 1, 1, 1e15)[[1]]
#makeCanonicalInitialization(list(list(theta=10, beta=2)), degree.sampled.vec, c(1,2,3), c(1,2), 1)[[1]]
#Njs<-rpois(20,10)
#names(Njs)<- seq(along.with=Njs)
#makeCanonicalInitialization(list(list(theta=10, beta=2, Njs=Njs)), degree.sampled.vec, c(1,2,3), c(1,2), 1)[[1]]









# Creates loose bound on maximal (log) beta values given theta and data:
BetaBound <- function(Js, Nj, theta) {
	# Njs assumed to be a table with degree (j) and counts (Njs)
	max.observed.degree<- max(Js)
	maximal.degree.count<- max(Nj)
	pop.size<- Nj %*% Js
	beta<- 1/ (pop.size * maximal.degree.count * max.observed.degree^theta )
	return(beta)
}
## Testing:
#BetaBound(sampled.degree.vector = rpois(200,2), theta = 1)














estimate.rds<- function (sampled.degree.vector, Sij, initial.values, arc=FALSE, control=generate.rds.control(), all.solutions=FALSE) {  	
	# Initializing:
	the.call<- sys.call()
	method <- control$method
	maxit <- control$maxit
	fnscale <- control$fnscale
	Nj.inflations<- control$Nj.inflations
	beta.inflations<- control$beta.inflations
	beta.scale <- control$beta.scale
	
	max.observed.degree<- max(sampled.degree.vector)
	N.j<- rep(0, max.observed.degree)
	# Look for the degrees in the data with non trivial estimates:
	Observed.Njs<-table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js<- as.numeric(names(Observed.Njs))
	maximal.degree.count<- max(Observed.Njs)
	final.result<- NA
	
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrap the likelihood function. 
	# Implements constraints on parameters by taking real valued (canonical) parameters and remapping them.
	likelihood.wrap<- function(par){
		# Initialize:
		likelihood.result<- -9999999
		
		beta<- inv.transform.beta(par[['canonical.beta']]) # assumes beta is non negative and given in log scale
		theta<- inv.transform.theta(canonical.theta = par[['canonical.theta']], control)
		
		N.j<- rep(0, max.observed.degree)
		get.j.index<- function(x) grep(paste("canonical.Nj.",x,"$",sep=""), names(par) ) # Get the indexes in par of the Njs.  
		observed.js.indexes<- 	sapply( Observed.js,  get.j.index)
		N.j[Observed.js]<- inv.transform.Nj(par[observed.js.indexes]) # fill non trivial Nj estimates.
		
		
		# Checking estimates are within allowed range:
		bad.Nj<- any( round(N.j[Observed.js],10) < round(Observed.Njs,10) ) # Estimated Nj smaller than observed
		bad.beta<- beta >= BetaBound(Js= Observed.js, Nj=N.j[Observed.js], theta=theta)  # Defined beta values
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
	
	
	# TODO: B) Automate initialization values for theta. (maybe using moments of intervals between samples).
	# TODO: A) Use better N_j initialization: IDW or inflated values.
	
	
	### Generate initial values:
	# Fills and formats initialization values as required by optim.wrap()
	# Assumes initial.values is a list of lists. Each containing a theta, beta and an Nj vector.	
	
	initial.value.list <- makeCanonicalInitialization(
			initial.values = initial.values, 
			sampled.degree.vector = sampled.degree.vector, 
			Nj.inflations = Nj.inflations, 
			beta.inflations = beta.inflations, 
			scale = beta.scale)
	
	# In case of convergence problems, try different optimization methods. In particular: "Nelder-Mead", "SANN",
	optim.wrap<- function(x) {
		try(optim(par=unlist(x), fn=likelihood.wrap, method=method , control=list(fnscale=fnscale, maxit=maxit)))
	}
	
	likelihood.optim<-lapply(initial.value.list,  optim.wrap )
	
	
	# Convert output from canonical form to original form:
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge")
		
		if(length(x)==5L){
			observed.js.indexes<- unlist(lapply(Observed.js, function(y) grep(paste("canonical.Nj.",y,"$",sep=""), names(x$par))))
			N.j[Observed.js]<- exp(x$par[observed.js.indexes]) # fill non trivial Nj estimates.		
			
			result<- list( 
					beta=inv.transform.beta(x$par[['canonical.beta']]), 
					theta=inv.transform.theta(x$par[['canonical.theta']], control),
					initial.values=initial.values,
					Nj=N.j,
					iterations=x$counts,
					likelihood.optimum=x$value, 
					call=the.call)			
		}		
		return(result)
	} 
	
	temp.result<- lapply(likelihood.optim, prepare.result)
	
	if(  any(sapply(temp.result, length) > 2 )   ){
		clean.temp.result<- as.list(temp.result[sapply(temp.result, length) > 2])
		if(!all.solutions){ 
			final.result<- temp.result[[  which.max(sapply(clean.temp.result, function(x) x$likelihood.optimum)) ]]
		}
		else{
			final.result<- clean.temp.result
		} 
	}
	else{
		stop('Estimation did not converge. Try differet intialization values or optimization method.')
	}
	
	return(final.result)							
}

##### Testing: 
#require(chords)
#data(simulation, package='chords')
#temp.data<- unlist(data3[1,7000:7500])
## Initialize only with thetas:
#initial.values<- list(list(theta=-1),list(theta=1))
### TODO: A)fix estimate.rds
#(rds.result<- estimate.rds(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), 
#					initial.values=initial.values, 
#					control=generate.rds.control(maxit=2000, 
#							method="Nelder-Mead", 
#							beta.inflations = exp(-seq(2,0, by=-0.02)), 
#							Nj.inflations = 1, 
#							beta.scale = 1, 
#							fnscale = -10)))
#plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')
#
#
#### Testing with Previous results:
#load(file="/Users/jonathanrosenblatt/Dropbox/Yakir/White/testingVer0.69.RData")
#.C("likelihood", 
#		sample=as.integer(data.degree), 
#		Sij=as.integer(as.matrix(data.Sjt)),
#		S=as.integer(compute.S(data.Sjt)),				  
#		beta=as.numeric(initial.values[[1]]$beta), 
#		theta=as.numeric(initial.values[[1]]$theta), 
#		Nj=as.numeric(initial.values[[1]]$Njs), 
#		observed_degrees=as.integer(rownames(data.Sjt)),
#		n=as.integer(length(data.degree)), 
#		N=as.integer(length(initial.values[[1]]$Njs)),
#		N_observed=as.integer(nrow(data.Sjt)),				  
#		arc=FALSE,
#		result=as.double(0))
#
#
#
#
#
#names(initial.values[[1]])[3]<-"Njs" 
#estimation4<- estimate.rds(sampled.degree.vector=data.degree, 
#		Sij = data.Sjt, 
#		initial.values=initial.values, 
#		control = generate.rds.control(
#				maxit=0, 
#				method="BFGS", 
#				beta.inflations = 1, 
#				Nj.inflations = 1, 
#				beta.scale = 1, 
#				fnscale = -1                                 
#		))

