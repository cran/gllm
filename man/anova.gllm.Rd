\name{anova.gllm}
\alias{anova.gllm}
\title{Summarize multiple results from gllm}
\description{
Compare likelihood ratio test statistics from multiple calls to \code{gllm}.
}
\usage{
\method{anova}{gllm}(object, \dots, test=c("Chisq","none")) 
}
\arguments{
  \item{object}{is an object output from gllm.}
  \item{\dots}{other objects from gllm.}
  \item{test}{evaluate LRTS for model, or nothing.}
}
\value{ 
  A list with components:
  \item{Model}{name of each object being compared}
  \item{Resid.df}{residual degrees of freedom for each model}
  \item{Deviance}{likelihood ratio test statistic for model versus saturated model}
  \item{Pr.Fit}{chi-square based P-value for model}
  \item{Test}{models compared in stepwise testing}
  \item{Df}{degrees of freedom of model comparson}
  \item{LRtest}{likelihood ratio test statistic comparing models}
  \item{Prob}{chi-square based P-value for LRTS}.
}
\author{ David L Duffy }
\keyword{ print }
