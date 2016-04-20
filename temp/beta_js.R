# 
# Author: johnros
###############################################################################



## TODO: B) Pad exported beta_js with NaNs or give name attribute to output.
## TODO: A) Initialize using final values of theta and beta.
## TODO: A) Lieklihood now expectes beta_js for ALL js (not only observed) update functions accordingly.


estimate.rds.free.betas<- function (sampled.degree.vector, Sij, method="BFGS", initial.values, arc=FALSE, control=generate.rds.control()) {
	##Important note: it is assumed initial.values has all parameter values and not just observable (thus estimable) ones! 
	
	# Initializing:
	max.observed.degree<- max(sampled.degree.vector)
	N.j<- rep(0, max.observed.degree)
	final.result<- NA
	# Look for the degrees in the datawith non trivial estimates:
	Observed.Njs<-table(sampled.degree.vector)[-1] # table of degree counts withuot zero!
	Observed.js<- as.numeric(names(Observed.Njs))
	maximal.degree.count<- max(Observed.Njs)
	the.call<- sys.call()
	
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrap the likelihood function. Implements constraints on parameters by taking real valued parameters and remapping them to the constrained space.
	likelihood.wrap<- function(par){
		# Initialize:
		likelihood.result<- -9999999 
		
		N.j<- rep(0, max.observed.degree)
		observed.Nj.indexes<- 	sapply(Observed.js, function(x) grep(paste("logNjs.",x,"$",sep=""), names(par) )   )
		N.j[Observed.js]<- exp(par[observed.Nj.indexes]) # fill non trivial Nj estimates.
		
		observed.beta_js.indexes<- 	sapply(Observed.js, function(x) grep(paste("canonical.betas.",x,"$",sep=""), names(par) )   )
		beta_js<- rep(0, max.observed.degree)
		beta_js[Observed.js]<- exp(par[observed.beta_js.indexes]) # fill non trivial Nj estimates.
				
		# Checking estimates are within allowed range:
		if(isTRUE(   any(round(N.j[Observed.js],10) < round(Observed.Njs,10))  ||   any(beta_js >= 1/(sum(N.j) * N.j))  ))   {
			return(likelihood.result)
		} 
		
		if(is.numeric(beta_js) &&  all(!is.infinite(beta_js)) ) {
			result<- .C("likelihood_beta_js", 
					sample=as.integer(sampled.degree.vector), 
					Sij=as.integer(Sij),
					S=as.integer(compute.S(Sij)),				  
					beta_js=as.numeric(beta_js), 
					Nj=as.numeric(N.j), 
					observed_degrees=as.integer(rownames(Sij)),
					n=as.integer(length(sampled.degree.vector)), 
					N=as.integer(length(N.j)),
					N_observed=as.integer(nrow(Sij)),				  
					arc=FALSE,
					result=as.double(0))			
		}
		if(!is.infinite(result$result)) likelihood.result<-result$result
		return(likelihood.result)
	}
	
	
	generate.initial.values<- function(initial.values){
		# maximal possible beta:
		log.beta_js.constraint<- log( 1/ (sum(sampled.degree.vector) * Observed.Njs )	)-0.01
		
		
		
		if(missing(initial.values)){
			wrap.initial.values<- function(...){
				c(
						canonical.betas=log.beta_js.constraint,						
						logNjs=log(Observed.Njs)+0.1
				)}
			returned.initial.values <- lapply(list(1), wrap.initial.values)			
		}
		
		else if(length(initial.values)==2L){
			# Note: Initial values are assumed to discard zeroes and include ALL parameters (not just observed)
			stopifnot(all(initial.values$beta_js >= 0))
			initial.log.betas_js<- log(initial.values$beta_js)			
			returned.initial.values<- list(c(
							canonical.betas=initial.log.betas_js,						
							logNjs=log(initial.values$Njs)+0.1 
							))			
		}
				
		
		else{
			stop("Unable to initialize optimization")
		}
		
				
		return(returned.initial.values)
	}
	
	
	
	
	
	
	
	## Initialize estimation:	 	
	initial.values.list<- generate.initial.values(initial.values)				
	
	# In case of convergence problems, try different optimization methods. In particular: "Nelder-Mead", "SANN", 
	optim.wrap<- function(x) try(optim(par=x, fn=likelihood.wrap, method=method , control=list(fnscale=-1, maxit=control$maxit)))
	
	likelihood.optim<-lapply(initial.values.list,  optim.wrap )
	
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge") 
		if(length(x) > 2){
			observed.js.indexes<- sapply(Observed.js, function(y) grep(paste("logNjs.",y,"$",sep=""), names(x$par)))
			N.j[Observed.js]<- exp(x$par[observed.js.indexes]) # fill non trivial Nj estimates.		
			
			result<- list( 
					beta_js= c(exp(x$par[-observed.js.indexes])),
					observed_js=Observed.js, 
					Nj=N.j,
					iterations=x$counts,
					likelihood.optimum=x$value, 
					call=the.call)			
		}		
		return(result)
	} 
	
	temp.result<- lapply(likelihood.optim, prepare.result)
	
	
	output.ok<- any(sapply(temp.result, length) > 2)	
	if(output.ok) {
		clean.temp.result<- temp.result[sapply(temp.result, length) > 2 ]
		final.result<- temp.result[[which.max(sapply(clean.temp.result, function(x) x$likelihood.optimum))]]
	}
	else{
		message('Estimation did not converge. Try differet intialization values or optimization method.')
	}
	
	return(final.result)							
}

##### Testing: 
#require(chords)
#data(simulation, package='chords')
#temp.data<- unlist(data3[1,7000:7500])
## Initialize only with thetas:
#(rds.result<- estimate.rds.free.betas(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), method="BFGS", 
#					control=generate.rds.control(maxit = 200)))
#str(rds.result)
#plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')
#x11()
#plot(rds.result$beta_js ~ sort(unique(temp.data))[-1])
#
#			 
#
#### Good values1: network.size = 500, network.density = 0.1, theta<- 3, beta<- 1e-10
#createDegreeCount <- function(network.size, network.density) {
#	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
#	return(table(unlist(neighbour.count[neighbour.count>0])))
#}
#(Njs<- createDegreeCount(network.size = 500, network.density = 0.1))
#theta<- 3
#beta<- 1e-10
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=7000))
#plot(degree.sampled.vec, type='h', main='Sample')
#
#
#
#
#
#### Manual creation:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
## Good values2: theta<- 3, beta<- 5e-10
#theta<- 1
#beta<- 5e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#plot(degree.sampled.vec, type='h', main='Sample')
#
#
#
#
#
## Testing with initialization from **true** values:
#true.js<- sort(unique(as.numeric(names(Njs))))
#(true.beta_js<- beta*true.js^theta)
#names(true.beta_js)<- true.js
#
#str(rds.result<- estimate.rds.free.betas(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), 
#				method='BFGS',				
#				initial.values = list(
#						beta_js=true.beta_js, 
#						Njs=Njs),  
#				control = generate.rds.control( maxit = 20)))
#{
#	plot(Njs, type='h', lwd=2.5)	
#	points(rds.result$Nj, type='h', col='orange', lwd=2)
#}
#{
#	plot(rds.result$beta_js~rds.result$observed_js)
#	lines( beta*rds.result$observed_js^theta  ~ rds.result$observed_js, lty=2)
#	lines(lowess(rds.result$beta_js~rds.result$observed_js), col='red')
#	lm.1<- lm(log(rds.result$beta_js)~log(rds.result$observed_js)-1)
#	lines(exp(predict(lm.1))~rds.result$observed_js , col='blue')
#	require(MASS)
#	lm.2<- rlm(log(rds.result$beta_js)~log(rds.result$observed_js)-1, method = "MM")
#	lines(exp(predict(lm.2))~rds.result$observed_js , col='green')
#}
#
#
#
#
#
#
#
#
#
### Testing without initialization:
#(rds.result<- estimate.rds.free.betas(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), 
#				method='BFGS',				
#				control = generate.rds.control( maxit = 200)))
#{
#	plot(rds.result$beta_js~rds.result$observed_js)
#	lines( beta*rds.result$observed_js^theta  ~ rds.result$observed_js, lty=2)
#	lines(lowess(rds.result$beta_js~rds.result$observed_js), col='red')
#	lm.1<- lm(log(rds.result$beta_js)~log(rds.result$observed_js)-1)
#	lines(exp(predict(lm.1))~rds.result$observed_js , col='blue')
#	require(MASS)
#	lm.2<- rlm(log(rds.result$beta_js)~log(rds.result$observed_js)-1, method = "MM")
#	lines(exp(predict(lm.2))~rds.result$observed_js , col='green')
#}
#{
#	plot(Njs, type='h', lwd=2.5)	
#	points(rds.result$Nj, type='h', col='orange', lwd=2)
#}





prepare.initialization<- function(free.beta.fit){
	zero.padded.Njs<- free.beta.fit$Nj
	zero.padded.js<- seq(along.with=zero.padded.Njs)
	zero.padded.beta_js<- free.beta.fit$beta * zero.padded.js^free.beta.fit$theta
	
	non.zero.indicators<- zero.padded.Njs>0
	Njs<- zero.padded.Njs[non.zero.indicators]
	names(Njs)<- zero.padded.js[non.zero.indicators]
	beta_js<- zero.padded.beta_js[non.zero.indicators]
	names(beta_js)<- zero.padded.js[non.zero.indicators]
	
	output<- list(beta_js=beta_js, Njs=Njs)
	return(output)
}
### Testing:
#require(chords)
#data(simulation, package='chords')
#temp.data<- unlist(data3[1,7000:7500])
#(rds.result<- estimate.rds(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), initial.values=list(theta=c(1,2,10)), 
#					control=generate.rds.control(maxit = 10), all.solutions = TRUE))
#prepare.initialization(rds.result[[1]])







estimate.rds.two.stage<- function(sampled.degree.vector, Sij, initial.values, arc=FALSE, control=generate.rds.control(), all.solutions=FALSE){
	method<- control$method
	first.estimates<- estimate.rds(sampled.degree.vector, Sij, method, initial.values, arc, control, all.solutions = TRUE)
	second.stage<- list()
	for(i in seq(along=first.estimates)){
		second.initialization<- prepare.initialization(first.estimates[[i]])
		second.stage[[i]]<- estimate.rds.free.betas(sampled.degree.vector, Sij, method, second.initialization, arc, control)
	}
	
	if(  any(sapply(second.stage, length) > 2 )   ){
		clean.second.stage<- second.stage[sapply(second.stage, length) > 2]
		if(!all.solutions){ 
			final.result<- second.stage[[  which.max(sapply(clean.second.stage, function(x) x$likelihood.optimum)) ]]
		}
		else{
			final.result<- clean.second.stage
		} 
	}
	else{
		message('Estimation did not converge. Try differet intialization values or optimization method.')
	}
	
	return(final.result)
}
## Testing:
# Good values1: network.size = 500, network.density = 0.1, theta<- 3, beta<- 1e-10
#createDegreeCount <- function(network.size, network.density) {
#	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
#	return(table(unlist(neighbour.count[neighbour.count>0])))
#}
#{
#	Njs<- createDegreeCount(network.size = 500, network.density = 0.1)
#	theta<- 3
#	beta<- 1e-10
#	tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1000))
#	plot(degree.sampled.vec, type='h', main='Sample')
#	
#	str(rds.result<- estimate.rds.two.stage(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), 
#					method='BFGS',				
#					initial.values = list(theta=c(theta), beta=c(beta), Njs=list(Njs)),  
#					control = generate.rds.control( maxit = 10), all.solutions = FALSE))
#}
#{
#	plot(rds.result$beta_js~rds.result$observed_js)
#	lines(lowess(rds.result$beta_js~rds.result$observed_js), col='red')
#	lines( beta*rds.result$observed_js^theta  ~ rds.result$observed_js, lty=2)
#	(lm.1<- lm(log(rds.result$beta_js)~log(rds.result$observed_js)))
#	lines(exp(predict(lm.1))~rds.result$observed_js , col='blue')
#	require(MASS)
#	(lm.2<- rlm(log(rds.result$beta_js)~log(rds.result$observed_js), method = "MM"))
#	lines(exp(predict(lm.2))~rds.result$observed_js , col='green')
#}
#plot(col='lightgrey', Njs, type='h')
#
#
## Manual creation:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
# Good values2: theta<- 3, beta<- 5e-10
#theta<- 1
#beta<- 5e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#plot(degree.sampled.vec, type='h', main='Sample')



