\name{ld2.model}
\alias{ld2.model}
\title{Write design and filter matrices for log-linear model of
  linkage disequilibrium between two codominant autosomal loci}
\description{Write design and filter matrices for log-linear model of
  linkage disequilibrium between two codominant autosomal loci.
}
\usage{
ld2.model(nall1, nall2, formula="~a1+a2+p1+p2+d") 
}
\arguments{
  \item{nall1}{is number of alleles at first codominant locus.}
  \item{nall2}{is number of alleles at first codominant locus.}
  \item{formula}{is character string listing terms to be 
     included in model, where a1 denotes allele frequencies for locus 1,
     p1 the deviation from Hardy-Weinberg expectations for locus1,
     and d the intragametic allelic association parameters.}
}
\value{ 
  A list with components:
  \item{Geno}{is a dummy contingency table showing the expected order.}
  \item{s}{is a vector of indices, one for each cell of the full (unobserved) 
     contingency table, representing the appropriate cell of \code{y}}
  \item{X}{is the design matrix.}
}

\keyword{category}
\keyword{models}
