\name{summary.gllm}
\alias{summary.gllm}
\title{Summarize results of gllm}
\description{
Summarizes contents of result of call to \code{gllm}.
}
\keywords{ print
}
\usage{
summary.gllm(object,...) 
}
\arguments{
  \item{object}{is the object output from gllm.}
  \item{...}{other arguments}
}
\values{ 
  A list with components:
  \item{nobs} {the number of cells in the observed table}
  \item{nfull} {the number of cells in the full table}
  \item{mean.cell} {the mean cell count in the observed table}
  \item{deviance} {the final model deviance (-2 log likelihood)}
  \item{model.df} {the model degrees of freedom}
  \item{coefficients} {the model parameter estimates, standard errors}
  \item{residuals} {Observed and fitted counts, plus Pearsonian residuals}
}
