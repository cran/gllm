\name{emgllmfitter}
\alias{emgllmfitter}
\title{ Generalized log-linear modelling by EM and iterative proportional
fitting}
\description{
Fits log-linear models for incomplete contingency tables,
via an EM approach.
}
\usage{
emgllmfitter(y,s,X,maxit,tol) 
}
\arguments{
  \item{y}{is the observed contingency table.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
  \item{X}{is the design matrix.}
  \item{maxit}{is the number of EM iterations.}
  \item{tol}{is the convergence criterion for the LR criterion.}
}
\value{ 
  A list with components:
  \item{y}{the observed table}
  \item{ji}{s, the scatter vector}
  \item{c}{the design matrix}
  \item{istop}{maximum EM iterations}
  \item{conv}{the convergence tolerance}
  \item{e}{expected counts for the full (unobserved) table}
  \item{ni}{nrow(X)}
  \item{nj}{length(y)}
  \item{nk}{ncol(X)-1}
  \item{f}{expected counts}
}
\details{ 
The call to Andreas Borg's C code that fits the model by EM/IPF.  The
algorithm follows the approach in Haber (1984).
}
\references{
Haber M (1984).  AS207: Fitting a general log-linear model.  \emph{Appl
        Statist} 33:358-362.
}
\keyword{category}
\keyword{models}

