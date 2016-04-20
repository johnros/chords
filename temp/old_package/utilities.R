# 
# Author: jonathanrosenblatt
###############################################################################


#' Utility functions for creating snowball matrix.
#' Assumes input is a table of degrees and degree counts.
#' @param data The sample
#' @return The sum of degrees for each rank 
#' @author Jonathan Rosnblatt
#' @export
ranks.sum<- function(data){
	data.table<- table(data)
	result<- as.numeric(rownames(data.table))*c(data.table)
	return(result)
}



#' Make the Sij matrix assuming no droouts from the snowball
#' @param data 
#' @return TBC
#' @author johnros
#' @export
make.Sij<- function(data){
	stopifnot(length(unique(data))>1L)
	result<- matrix(0, ncol=length(data), nrow=max(data))
	for (i in 1:(length(data)-1)){ result[data[i],i+1]<- result[data[i],i+1]+1 }  
	result<- t(apply(result, 1 , cumsum))
	non.empty.ind<- apply(result, 1, sum)!=0
	.row.names<- as.character(seq(along.with=non.empty.ind)[non.empty.ind])
	result<- result[non.empty.ind,]
	rownames(result)<- .row.names
	return(result)
}






#' Compares the output of the optimizaton for different initialization values.
#' @param Sij 
#' @return Used for creating the Snowball size given an Sij matrix
#' @author johnros
#' @export
compute.S<- function(Sij){
	result<- apply(Sij, 2, sum)
	non.null.ind<- cumsum(result!=0)
	result[non.null.ind==0]<- 1    
	return(as.integer(result))
}







createDegreeCount <- function(network.size, network.density) {
	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
	return(table(unlist(neighbour.count[neighbour.count>0])))
}
## Testing: 
#(Njs<- createDegreeCount(network.size = 1e5, network.density = 0.01))








generate.sample <- function(theta, Njs, beta, sample.length, double.sampling.warning.thresh=0.05) {
	# Initializing:
	js<- as.numeric(names(Njs))
	N<- sum(Njs)
	maximal.beta<- 1/ (N * max(Njs) * max(js)^theta )
	double.sampling.warnings<- FALSE
	
	
	# Input validations:
	error.messge<- paste("beta value too large. Errors might occur. \nFor the given degree frequencies, consider keeping it under ", signif(maximal.beta), sep="" )
	if(beta > maximal.beta) message(error.messge) # Note: larger beta favour sampling non 0 degrees.
	stopifnot(min(js)>0)
	
	
	
	# Preparing to sample:
	Uik<- rep(0L, length.out=max(js))
	names(Uik)<- seq(along.with=Uik)
	snowball<- 1L
	degree.sampled.vec<- rep(NA, sample.length)
	# Start sampling:
	for(i in seq(1, sample.length)){	
		pi.function<- function(k) beta * k^theta * snowball * (Njs[paste(k)] - Uik[k] )
		n.pi.function<- function(k) 1-pi.function(k)
		# Compute the probabilities of sampling degree k=0:
		no.sample.probability<- prod(unlist(lapply(js, n.pi.function )))
		sample.someone<- as.logical(rbinom(1, 1, prob = 1 - no.sample.probability ))		
		if(sample.someone){
			# Compute the probabilities of sampling degree k>0:
			probs.function<- function(k){
				nominator<-  pi.function(k) * prod(unlist(lapply(js[js!=k], n.pi.function)))  
				# Note: the true probability is not needed due to the normalization in sample(). A constant is used insted.
				# denominator<-  1 - prod(unlist(lapply(js, n.pi.function)))
				denominator<-  1 			
				return(nominator/denominator)		
			}
			degree.probabilities<- unlist(lapply(js, probs.function))			
			
			if( sum(degree.probabilities)+ no.sample.probability <  1- double.sampling.warning.thresh) double.sampling.warnings<- TRUE				
			degree.sampled<- sample(x=js, prob=degree.probabilities, size=1)      
			Uik[degree.sampled]<- Uik[degree.sampled]+1L
			snowball<- sum(Uik)
			degree.sampled.vec[i]<- as.integer(degree.sampled)    
		}  
		else{ 
			degree.sampled.vec[i]<- 0L
		}
	}
	
	# Exiting:
	if(double.sampling.warnings) message("Double sampling probabilities non-negligeable. Consider lower beta values.")
	return(degree.sampled.vec)
}
##Testing:
#Njs<- c(100,500,500,200)
#names(Njs)<- c("1","50","100","1000")
#Njs<- as.table(Njs)
#theta<- 1.1
#beta<- 3e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#(.Nj.table <- table(degree.sampled.vec))




var.theta<- function(sampled.degree.vector, Sij=make.Sij(sampled.degree.vector), Njs, theta, beta){
	I<- compute.S(Sij)
	N<- sum(Njs*seq(along.with=Njs))
	Nj.uniques<- sort(unique(sampled.degree.vector))
	total.sampling<- length(sampled.degree.vector)
	make.value<- Vectorize(function(k,t){
				if(!(k %in% rownames(Sij))) return(0)
				else return(
							log(k)^2* k^theta* I[t]/N *(Njs[k]/N - Sij[as.character(k),t]/N)
					)
			})
	information<- beta * N * sum(outer(Nj.uniques, seq(1,total.sampling), FUN="make.value"))
	return( (1/information)/N )
}
## Example:
#var.theta(degree.sampled.vec, make.Sij(degree.sampled.vec), rds.result$Nj, rds.result$theta, rds.result$beta )



