\name{boot.table}
\alias{boot.table}
\title{Produce one bootstrap replicate of a vector of counts}
\description{
Given a vector of counts from a contingency table, produce a bootstrap 
replicate.  Sampling zeroes are replaced by 0.5.
}
\keywords{ manip
}
\usage{
boot.table(y) 
}
\arguments{
  \item{y}{is the observed contingency table.}
 }
}
\values{ 
  A vector of counts with the same total.
}
\examples{
boot.table(c(1,3,4,2))
## 0.5 2 5 3
}
