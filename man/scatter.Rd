\name{scatter}
\alias{scatter}
\title{Create a filter matrix from a summary array of indices}
\description{
Create a filter matrix that multiplying the vector of counts from a 
complete contingency table, gives a collapsed contingency table.
}
\keywords{ manip
}
\usage{
scatter(y,s) 
}
\arguments{
  \item{y}{is the observed contingency table. Provides a target length only.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
 }
}
\values{ 
  \item{S} {A matrix of orthogonal columns of 1s and 0s}
}
\examples{
y<-double(3)
z<-1:5
z %*% scatter(y,c(1,1,2,3,3)
## 1+2, 3, 4+5
}
