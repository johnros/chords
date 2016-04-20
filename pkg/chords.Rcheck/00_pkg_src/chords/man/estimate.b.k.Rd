\name{estimate.b.k}
\alias{estimate.b.k}

\title{
Estimate population size from RDS sample
}
\description{
Performs maximum likelihood estimation of a population size from a respondent driven sample.
The underlying idea is that the sample spreads like an epidemic in the target population as described in the reference.
}
\usage{
  estimate.b.k(rds.object, 
  const=1,
  impute.Nks=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds.object}{
  A object of class \code{rds-object} as constructed by \code{\link{initializeRdsObject}}
}
\item{const}{A scaler value for numerical stability. Should not be altered.}
\item{impute.Nks}{Should inestimable values be imputed?}
}

\details{
See the reference for a description of the likelihood problem solved.
Optimization is performed by noting that likelihood is coordinate-wise convex, thus amounts to a series of line-searches.

If \code{impute.Nk} is enabled (default), a simpe heuristic is used to estimate degree frequencies which are inestimable using maximum likelihood. This essentially consists of scaling the observed degree by the average observed-to-estimated factor over all estimable frequencies.
}

\note{
This function replaces the \code{estimate.rds} function from versions prior to 0.8 which is now  deprecated.
Differences between the functions include:
1- Likelihood is now computed in continuous time. The sample data is thus assumed to have an interviewDt column.
2- Optimization now exploits the coordinate-wise convexity of the likelihood function for efficient solving.
}

\value{
A list with the following components:
\item{call}{The function call.}
\item{Nk.estimates}{The estimated degree frequencies.}
\item{log.bk.estimates}{The estimated sampling rates for each degree.}
}


\references{
  [1] Berchenko, Yakir, Jonathan Rosenblatt, and Simon D. W. Frost. "Modeling and Analyzing Respondent Driven Sampling as a Counting Process."" arXiv:1304.3505, April 11, 2013. HTTP://arXiv.org/abs/1304.3505.
  }

\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}

\seealso{
  \code{\link{initializeRdsObject}}, \code{\link{makeRdsSample}}, \code{\link{getTheta}}
}
\examples{
# Import data in RDS file format:
data(brazil)
# Initialize RDS object:
rds.object2<- initializeRdsObject(brazil)
# Estimate:
rds.object2$estimates <- estimate.b.k(rds.object = rds.object2 )

# View estimates:
plot(rds.object2$estimates$Nk.estimates, type='h')
# Population size estimate:
sum(rds.object2$estimates$Nk.estimates)
plot(rds.object2$estimates$log.bk.estimates, type='h')


## Recover theta assuming b.k=b_0*k^theta
getTheta(rds.object2)

# How many degrees were imputed?:
table(rds.object2$estimates$convergence)
}
