\name{boot.table}
\alias{boot.table}
\title{Produce one bootstrap replicate of a vector of counts}
\description{
Given a vector of counts from a contingency table, produce a bootstrap 
replicate.  Sampling zeroes are replaced by 0.5.
}
\usage{
boot.table(y,strata=NULL) 
}
\arguments{
  \item{y}{is the observed contingency table.}
  \item{strata}{is a vector defining the strata for a stratified bootstrap.}
}
\value{ 
  A vector of counts with the same total.
}
\examples{
boot.table(c(1,3,4,2))
## 0.5 2.0 5.0 3.0
boot.table(c(1,3,4,2),c(1,2,1,2))
## 2 1 3 4
}
\keyword{category}
\keyword{manip}
