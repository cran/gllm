\name{scoregllm}
\alias{scoregllm}
\title{ Generalized log-linear modelling via Fisher scoring }
\description{
Fits log-linear models for incomplete contingency tables,
including some latent class models, via Fisher scoring approaches.
}
\usage{
scoregllm(y,s,X,m,tol=1e-5)
}
\arguments{
  \item{y}{is the observed contingency table.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
  \item{X}{is the design matrix or a formula.}
  \item{m}{is a vector of starting values for the full (unobserved)
     contingency table.}
  \item{tol}{is the convergence criterion for the LR criterion.}
}
\values{ 
  A list with components:
  \item{iter} {the number of scoring iterations until convergence}
  \item{deviance} {the final model deviance (-2 log likelihood)}
  \item{df} {the model degrees of freedom}
  \item{coefficients} {the model parameter estimates}
  \item{se} {the standard errors for the model parameter estimates}
  \item{V} {the variance-covariance matrix for the model parameter estimates}
  \item{observed.values} {the observed counts in \code{y}}
  \item{fitted.values} {the expected counts under the fitted model}
  \item{residuals} {Pearsonian residuals under the fitted model}
  \item{full.table} {the expected counts for the full (unobserved) table.}
}
\details{ 
The generalized log-linear model allows for modelling of incomplete
contingency tables, that is tables where one or more dimensions have
been collapsed over.  These include situations where imprecise
measures have been calibrated using a "perfect" gold standard, and the
true association between imperfectly measured variables is to be
estimated; where data is missing for a subsample of the population;
latent variable models where latent variables are "errorless" functions
of observed variables - eg ML gene frequency estimation from counts of
observed phenotypes; specialised measurement models eg where observed
counts are mixtures due to perfect measures and error prone measures;
standard latent class analysis; symmetry and quasi-symmetry models for
square tables.
  
The general framework underlying these models is summarised by Espeland
(1986), and Espeland & Hui (1987), and is originally due to Thompson &
Baker (1981).  An observed contingency table \eqn{y}, which will be
treated as a vector, is modelled as arising from an underlying complete
table \eqn{z}, where observed count \eqn{y_j} is the sum of a number of
elements of \eqn{z}, such that each \eqn{z_i} contributes to no more
than one \eqn{y_j}. Therefore one can write \eqn{y=F'z}, where \eqn{F}
is made up of orthogonal columns of ones and zeros.

We then specify a loglinear model for \eqn{z}, so that
\eqn{log(E(z))=X'b}, where \eqn{X} is a design matrix, and \eqn{b} a
vector of loglinear parameters.  The loglinear model for \eqn{z} and
thus \eqn{y}, can be fitted by a Fisher scoring approach, presented in
Espeland (1986).

The \code{gllm} function is actually a simple wrapper for \code{scoregllm()}.
}
\references{
Espeland MA (1986).  A general class of models for discrete multivariate
        data. \emph{Commun. Statist.-Simula} 15:405-424.

Espeland MA, Hui SL (1987).  A general approach to analyzing epidemiologic
        data that contains misclassification errors.  \emph{Biometrics}
        43:1001-1012.

Thompson R, Baker RJ (1981).  Composite link functions in generalized
        linear models.  \emph{Appl Statist} 30: 125-131.
}

\examples{ 
#
# An example of model fitting: gametic association between two diallelic loci
# Data matrix
#
y<-c( 187,386,156,
      352,310,20,
      136,0  ,0)
#
# Scatter matrix
#
s<-  c( 1, 2, 2, 3,
        4, 5, 5, 6,
        4, 5, 5, 6,
        7, 8, 8, 9)
#
# Design matrix
#
X<-  matrix(c( 1,0,0,0,0,0,1,
               1,0,1,0,0,0,0,
               1,0,1,0,0,0,0,
               1,0,2,0,1,0,0,
               1,1,0,0,0,0,0,
               1,1,1,0,0,1,0,
               1,1,1,0,0,0,1,
               1,1,2,0,1,1,1,
               1,1,0,0,0,0,0,
               1,1,1,0,0,0,1,
               1,1,1,0,0,1,0,
               1,1,2,0,1,1,1,
               1,2,0,1,0,0,0,
               1,2,1,1,0,1,1,
               1,2,1,1,0,1,1,
               1,2,2,1,1,2,2), byrow=TRUE, ncol=7) 

colnames(X)<-c("Intercept", "A", "B", "P1", "P2", "Delta", "Epsilon")
res<-scoregllm(y,s,X[,c(1:6)], 
               c(255,176,176,121,164,37,113,25,164,113,37,25,90,20,20,5))
summary(res)
#

}
\keyword{ models }
\keyword{ category }
