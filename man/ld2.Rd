\name{ld2}
\alias{ld2}
\alias{print.ld2}
\title{Estimate linkage disequilibrium between two codominant autosomal loci}
\description{
Fits a log-linear model for allelic association between two codominant
autosomal loci.  Measures of LD are odds ratios.
}
\usage{
ld2(locus1, locus2) 
}
\arguments{
  \item{locus1}{is a character vector containing the genotypes at the first
     locus, or a RxC contingency table of genotype counts.}
  \item{locus2}{is a character vector containing the genotypes at the second
     locus.}
}
\values{ 
  \item{m0} {base model}
  \item{m1} {estimating LD coefficient(s) assuming HWE}
  \item{m2} {testing HWE at locus 1}
  \item{m3} {testing HWE at locus 2}
  \item{m4} {estimating LD and HWD coefficient(s)}
}
\examples{
MNS<-matrix(c(91,32,5,147,78,17,85,75,7), nr=3)
colnames(MNS)<-c("S/S","S/s","s/s")
rownames(MNS)<-c("M/M","M/N","N/N")
class(MNS)<-"table"
print(MNS)
res<-ld2(MNS)
print(res)
}
\keyword{category}
\keyword{models}
