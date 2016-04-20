\name{thetaSmoothingNks}
\alias{thetaSmoothingNks}

\title{
  Smooth estimates degree frequencies.
}
\description{
  Smoothes estimated \eqn{N_k} by assuming that \eqn{beta_k=beta_0 * k^theta}.
}
\usage{
  thetaSmoothingNks(rds.object,...)
}

\arguments{
  \item{rds.object}{A \code{rds-object} class object as returned by \code{\link{estimate.b.k}}}
  \item{...}{Firther arguments passed to \code{\link{getTheta}}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A numeric vector of smoothed \eqn{N_k} values.}

\author{
Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}

\seealso{
\code{\link{estimate.b.k}}; \code{\link{getTheta}}
}

\examples{
# See estimate.b.k()
}
