\name{emgllm}
\alias{emgllm}
\title{ Generalized log-linear modelling by EM and iterative proportional
fitting}
\description{
Fits log-linear models for incomplete contingency tables,
including some latent class models, via an EM approach.
}
\usage{
emgllm(y,s,X,maxit=1000,tol=0.00001) 
}
\arguments{
  \item{y}{is the observed contingency table.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
  \item{X}{is the design matrix, or a formula.}
  \item{maxit}{is the number of EM iterations.}
  \item{tol}{is the convergence criterion for the LR criterion.}
}
\values{ 
  A list with components:
  \item{deviance} {the final model deviance (-2 log likelihood)}
  \item{observed.values} {the observed counts in \code{y}}
  \item{fitted.values} {the expected counts under the fitted model}
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
thus \eqn{y}, can be fitted via an iterative proportional fitting algorithm
for \eqn{b} and \eqn{z}, with an EM fitting for \eqn{y}, \eqn{z} and
\eqn{b} (Haber 1984).  

The \code{emgllm} function is a wrapper for the slightly modified AS207
Fortran code from Haber (1984).
}
\references{
Espeland MA (1986).  A general class of models for discrete multivariate
        data. \emph{Commun. Statist.-Simula} 15:405-424.

Espeland MA, Hui SL (1987).  A general approach to analyzing epidemiologic
        data that contains misclassification errors.  \emph{Biometrics}
        43:1001-1012.

Haber M (1984).  AS207: Fitting a general log-linear model.  \emph{Appl
        Statist} 33:358-362.

Thompson R, Baker RJ (1981).  Composite link functions in generalized
        linear models.  \emph{Appl Stat} 30: 125-131.
}

\examples{ 
#
# latent class analysis: two latent classes
#
# Data matrix 2x2x2x2x2 table of responses to five binary items
#
y<-c( 3,   6,   2,   11,   1,   1,   3,    4,
      1,   8,   0,   16,   0,   3,   2,   15,
     10,  29,  14,   81,   3,  28,  15,   80,
     16,  56,  21,  173,  11,  61,  28,  298)
#
# Scatter matrix: full table is 2x2x2x2x2x2
#
s<-  c(1:32,1:32)
#
# Design matrix: x is the latent variable (2 levels), 
# a-e are the observed variables
#
i<-rep(1,64)
x<-as.integer(gl(2,32,64))-1
a<-as.integer(gl(2,16,64))-1
b<-as.integer(gl(2,8 ,64))-1
c<-as.integer(gl(2,4 ,64))-1
d<-as.integer(gl(2,2 ,64))-1
e<-as.integer(gl(2,1 ,64))-1
X<-cbind(i,x,a,b,c,d,e,x*cbind(a,b,c,d,e))
colnames(X)<-c("Int","X","A","B","C","D","E","AX","BX","CX","DX","EX")
res<-emgllm(y,s,X, tol=0.01)
res
#
# Obtain standard errors for parameter estimates
#
summary(scoregllm(y,s,X,as.array(res$full.table)))
}
\keyword{category}
\keyword{models}
