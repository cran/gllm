\name{boot.gllm}
\alias{boot.gllm}
\title{ Bootstrap for generalized log-linear modelling }
\description{
Fits log-linear models for incomplete contingency tables,
including some latent class models, via EM and Fisher scoring approaches.
Performs a bootstrap for the sampling distribution of the 
full unobserved table.
}
\usage{
boot.gllm(y,s,X,method="hybrid",em.maxit=1,tol=0.00001,strata=NULL,R=200)
}
\arguments{
  \item{y}{is the observed contingency table.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
  \item{X}{is the design matrix, or a formula.}
  \item{method}{chooses the EM, Fisher scoring or a hybrid (EM then scoring)
     method for fitting the model.}
  \item{em.maxit}{is the number of EM iterations.}
  \item{tol}{is the convergence criterion for the LR criterion.}
  \item{strata}{is a vector identifying the sampling strata.}
  \item{R}{is the number of bootstrap replicates.}
}
\values{ 
  A matrix \eqn{R+1} by \code{ncol(X)} containing the initial estimate
  of the full (unobserved) contingency table, and the \eqn{R} bootstrap
  replicates of the full table.
}
\details{ 
The generalized log-linear model allows for modelling of incomplete
contingency tables, that is tables where one or more dimensions have
been collapsed over.  See \code{gllm} for details.

Often, functions of the full unobserved table are the main focus of the
analysis.  For example, in a double sampling design where there is a
gold standard measure for one part of the data set and only an
unreliable measure for another part, the expected value of the gold
standard in the entire dataset is the outcome of interest.  The standard
error of this statistic may be a complex function of the observed counts
and model parameters.

Bootstrapping is one way to estimate such standard errors from a complex
sampling design.  The bootstrap sampling may be stratified if the design
implies this, e.g. product-multinomial.

}
\references{
Hochberg Y (1977).  On the use of double sampling schemes in analyzing
        categorical data with misclassification errors.  \emph{J Am Statist
        Ass} 72:914-921. 
}

\examples{ 
#
# Fit Hochberg 1977 double sampling data
# 2x2 table of imprecise measures and 2x2x2x2 reliability data
#
# 2x2 table of imprecise measures
#
y1 <-c(1196, 13562,
       7151, 58175)
a2<- 2-as.integer(gl(2,1,4))
b2<- 2-as.integer(gl(2,2,4))
set1<-data.frame(y1,a2,b2)
#
# 2x2x2x2 reliability data
#
y2<-c(17, 3,   10, 258,
       3, 4,    4,  25,
      16, 3,   25, 197,
     100, 13, 107, 1014)

a <- 2-as.integer(gl(2,1,16))
a2<- 2-as.integer(gl(2,2,16))
b <- 2-as.integer(gl(2,4,16))
b2<- 2-as.integer(gl(2,8,16))

set2<-data.frame(y2,a,a2,b,b2)
#
# Combined analysis
#
y<-c(y1,y2)
#
# Map observed table onto underlying 2x2x2x2x2 table
#
s <-c(1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4,
      5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
#
# Model combining the tables is A*A2*B*B2 + L (dummy study variable)
#
a <- 2-as.integer(gl(2,1,32))
a2<- 2-as.integer(gl(2,2,32))
b <- 2-as.integer(gl(2,4,32))
b2<- 2-as.integer(gl(2,8,32))
l <- 2-as.integer(gl(2,16,32))

X <- model.matrix( ~ a*a2*b*b2+l)

#
# Table 1 using unreliable measure
#
res1<-glm(y1 ~ a2*b2, family=poisson(),data=set1)
print(summary(res1))
#
# Table 2 using reliable measure
#
res2a<-glm(y2 ~ a*b, family=poisson(),data=set2)
print(summary(res2a))
#
# Table 2 demonstrating complex relationship between gold standard and
# unreliable measure
#
res2b<-glm(y2 ~ a*a2*b*b2, family=poisson(),data=set2)
print(summary(res2b))
#
# Combined analysis
#
require(gllm)
res12<-gllm(y,s,X)
print(summary.gllm(res12))
#
# Bootstrap the collapsed table to get estimated OR for reliable measures
#
# a and b are binary vectors the length of the *full* table
# and define the variables for which the odds ratio is to be
# estimated, here the reliable measure of injury and seatbelt usage.
#
boot.hochberg <- function (y,s,X,nrep,a,b) {
  z<-boot.gllm(y,s,X,R=nrep)
  boot.tab<-cbind(apply(z[,a & b],1,sum),
                  apply(z[,!a & b],1,sum),
                  apply(z[,a & !b],1,sum),
                  apply(z[,!a & !b],1,sum))
  oddsr<-boot.tab[,1]*boot.tab[,4]/boot.tab[,2]/boot.tab[,3] 
  hochberg.tab<-data.frame( c("yes","yes","no","no"),
                            c("yes","no","yes","no"),
                            boot.tab[1,],
                            apply(boot.tab[2:(1+nrep),],2,sd))
  colnames(hochberg.tab)<-c("Precise Injury","Precise seatbelt usage",
                            "Estimated Count","Bootstrap S.E.")
  print(hochberg.tab)
  cat("\nEstimated OR=",oddsr[1],"\n")
  cat("        Bias=",oddsr[1]-mean(oddsr[2:(1+nrep)]),"\n")
  cat("Bootstrap SE=",sd(oddsr[2:(1+nrep)]),"\n\nQuantiles\n\n")
  quantile(oddsr[2:(1+nrep)],c(0.025,0.50,0.975))
}
boot.hochberg(y,s,X,nrep=20,a,b)
}
\keyword{category}
\keyword{models}
