\name{as207}
\alias{as207}
\title{ Generalized log-linear modelling by EM and iterative proportional
fitting}
\description{
Fits log-linear models for incomplete contingency tables,
via an EM approach.
}
\keywords{ models regression multivariate
}
\usage{
as207(y,s,X,maxit,tol) 
}
\arguments{
  \item{y}{is the observed contingency table.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
  \item{X}{is the design matrix.}
  \item{maxit}{is the number of EM iterations.}
  \item{tol}{is the convergence criterion for the LR criterion.}
 }
}
\values{ 
  A list with components:
  \item{istop} {maximum EM iterations}
  \item{ni} {nrow(X)}
  \item{nid} {nrow(X)}
  \item{nj} {length(y)}
  \item{nk} {ncol(X)-1}
  \item{nkp}{ncol(X)} 
  \item{ji} {s, the scatter vector}
  \item{y} {the observed table}
  \item{c} {the design matrix}
  \item{conv} {the convergence tolerance}
  \item{w} {work matrix}
  \item{v} {work matrix}
  \item{e} {expected counts for the full (unobserved) table}
  \item{f} {expected counts}
  \item{cspr} {Pearsonian \eqn{Chi^2}}
  \item{cslr} {Likelihood ratio \eqn{Chi^2}}
  \item{ifault} {error indicator}
}
\details{ 
The call to the slightly modified AS207 Fortran code from Haber (1984).
}
\references {
Haber M (1984).  AS207: Fitting a general log-linear model.  \emph{Appl
        Statist} 33:358-362.
}

