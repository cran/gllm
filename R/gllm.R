#
# Generalised log-linear modelling via EM and Fisher scoring
#
# Library setup
#
.First.lib <- function(lib,pkg) {
  library.dynam("gllm",pkg,lib)
  cat("This is gllm 0.20\n")
}
#
# EM IPF algorithm of Haber AS207
#
emgllm <- function(y,s,X,maxit=1000,tol=0.00001) {
  X<-cbind(X,double(nrow(X)))
  em<-as207(y,s,X,maxit,tol) 
  deviance<-em$cslr
  observed.values<-em$y
  fitted.values<-em$f
  full.table<-em$e
  return(deviance,
         observed.values,fitted.values,
         full.table)
}
as207 <- function(y,s,X,maxit,tol)
  .Fortran("gllm",
            istop = as.integer(maxit),
            ni = as.integer(nrow(X)),
            nid = as.integer(nrow(X)),
            nj = as.integer(length(y)),
            nk = as.integer(ncol(X)-1),
            nkp=  as.integer(ncol(X)),
            ji = as.integer(s),
            y = as.double(y),
            c = X,
            conv = as.double(tol),
            w = matrix(rep(0.0,4*nrow(X)),nrow=4),
            v = matrix(rep(0.0,2*(ncol(X))),nrow=2),
            e = double(nrow(X)),
            f = double(length(y)),
            cspr = double(1),
            cslr = double(1),
            ifault=as.integer(0))
#
# Convert scatter vector used by AS207 to matrix for scoring method
#
scatter<- function(y,s) {
  S<-matrix(rep(0,length(y)*length(s)),nrow=length(y))
  for(i in 1:length(s)) if (s[i]<=length(y)) S[s[i],i]<-1
  return(t(S))
}
#
# Espeland's scoring algorithm
# y=observed table 
# s=scatter vector 
# X=unaugmented design matrix
# m=starting values for complete table
#
scoregllm<-function(y,s,X,m,tol=1e-5) {
#
# Initialize
#
  eps<-0.00001
  X<-t(X)
  S<-scatter(y,s)
  z<-as.vector(S %*% solve(t(S) %*% S, tol=1e-10) %*% y)
  iter<- 0
  olddev<- -1
  deviance<- 0
  b<-solve(t(X),log(m),tol=1e-10)
#
# Main loop
#
  while(abs(olddev-deviance)>tol) {
    iter<-iter+1
    olddev<-deviance
    if (iter>1) m<- as.vector(exp(t(X) %*% b))
    P<- S %*% solve(t(S) %*% diag(m) %*% S, tol=1e-10) %*% t(S) %*% diag(m) 
    A<- P %*% t(X)
    V<- t(A) %*% diag(m) %*% A
    V<- qr.solve(V,tol=1e-10)
    b<- b - V %*% t(A) %*% (m - z)
    f<- as.vector(t(S) %*% m)
    use<-(y>eps & f>eps)
    deviance<- 2.0*sum(y[use]*log(y[use]/f[use]))
  }
  observed.values<-y
  fitted.values<-f
  residuals<- f
  residuals[f>0]<-(y[f>0]-f[f>0])/sqrt(f[f>0])
  full.table<-m
  coefficients<-as.vector(b)
  names(coefficients)<-rownames(X)
  bl<-(rownames(X)=="")
  names(coefficients)[bl]<-paste("beta",1:nrow(X),sep="")[bl]
  se<-sqrt(diag(V))
  df<-length(y)-qr(X)$rank
  res<-list(iter=iter,deviance=deviance,df=df,
            coefficients=coefficients,se=se,V=V,
            observed.values=observed.values,
            fitted.values=fitted.values,
            residuals=residuals,
            full.table=full.table)
  class(res) <- "gllm"
  res
}
#
# Global front end to call either/both routines
#
gllm <- function(y,s,X,method="hybrid",em.maxit=1,tol=0.00001) {
  if (method=="hybrid" || method=="scoring") {
    scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=em.maxit,tol=tol)$full.table)) 
  }else{
    scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=10000,tol=tol)$full.table)) 
  }
}
#
# Summary of results from gllm
#
summary.gllm <- function(object, ...) {
  tab.coef<-data.frame(object$coefficients, 
                       object$se, 
                       exp(object$coefficients),
                       exp(object$coefficients-1.96*object$se),
                       exp(object$coefficients+1.96*object$se),
                       row.names=names(object$coefficients))
  colnames(tab.coef)<-c("Estimate","S.E.","exp(Estimate)",  
                        "Lower 95% CL","Upper 95% CL")
  tab.fitted<-data.frame(object$observed.values, 
                         object$fitted.values, 
                         object$residuals)
  colnames(tab.fitted)<-c("Observed Count","Predicted","Residual")  
  summary<- list()
  summary$nobs<-length(object$observed.values)
  summary$nfull<-length(object$full.table)
  summary$mean.cell<-mean(object$observed.values)
  summary$deviance<-object$deviance 
  summary$model.df<-object$df
  summary$coefficients<-tab.coef
  summary$residuals<-tab.fitted
  class(summary) <- "summary.gllm"
  summary
}
#
# Print summary
#
print.summary.gllm <- function(x, digits=NULL, show.residuals = FALSE, ...) {
  if (is.null(digits))
    digits <- options()$digits
  else options(digits=digits)

  cat("\nNo. cells in observed table: ", x$nobs, "\n",sep="")
  cat("No. cells in complete table: ", x$nfull, "\n",sep="")
  cat("    Mean observed cell size: ", x$mean.cell, "\n",sep="")
  cat("        Model Deviance (df): ", formatC(x$deviance,digits=2,format="f"),
                                        " (", x$model.df, ")\n\n",sep="")

  print(x$coefficients, digits=digits)
  if (show.residuals) {
    cat("\n")
    print(x$residuals, digits=digits)
  }
}
#
# Bootstrap a contingency table.  Sampling zeroes augmented by 1/2.
#
boot.table <- function(y,strata=NULL) {  
  ynew<-rep(0.5, length(y))
  if (is.null(strata)) {
    tab.ynew<-table(sample(rep(1:length(y),y),replace=TRUE))
    ynew[as.integer(names(tab.ynew))]<-tab.ynew
  }else{
    s<-as.integer(strata)
    for(i in unique(s)) {
      idx<-s==i
      tab.ynew<-table(sample(rep((1:length(y))[idx],y[idx]),replace=TRUE))
      ynew[as.integer(names(tab.ynew))]<-tab.ynew
    }
  }
  ynew
}
#
# Bootstrap a GLLM returning the full fitted tables
#
boot.gllm <- function(y,s,X,method="hybrid",em.maxit=1,tol=0.00001,
                      strata=NULL,R=200) {
  if (method=="hybrid" || method=="scoring") {
    f0<-scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=em.maxit,tol=tol)$full.table))$full.table
  }else{
    f0<-scoregllm(y,s,X,as.array(emgllm(y,s,X,maxit=10000,tol=tol)$full.table))$full.table
  }
  result<-as.matrix(t(f0))
  cat("It",0,":",y,"\n")
  for(i in 1:R) {
    ynew<-boot.table(y,strata=strata)
    cat("It",i,":",ynew,"\n")
    result<-rbind(result, t(scoregllm(ynew,s,X,as.array(emgllm(ynew,s,X,
                            maxit=em.maxit,tol=tol)$full.table))$full.table))
  }
  result
}
#
# After anova.multinom
#
anova.gllm <- function(object, ..., test = c("Chisq", "none"))
{
  modelname<-unlist(strsplit(deparse(match.call()),"[(),]")) 
  modelname<-gsub("[() ]","",modelname)
  modelname<-gsub("object=","",modelname)
  modelname<-modelname[-c(1,grep("=",modelname))]
  print(modelname)
  test <- match.arg(test)
  dots <- list(...)
  if(length(dots) == 0)
    stop("anova is not implemented for a single gllm object")
  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$df)
  s <- order(dflis, decreasing=TRUE)
  mlist <- mlist[s]
  if(any(!sapply(mlist, inherits, "gllm"))) {
    stop("not all objects are of class `gllm'")
  }
  mds <- modelname[s]
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) x$deviance)
  tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model = mds, Resid.df = dfs,
                    Deviance = lls, Pr.Fit = 1-pchisq(lls,dfs), 
                    Test = tss, Df = df, LRtest = x2,
                    Prob = pr)
  names(out) <- c("Model", "Resid. df", "Resid. Dev", "Pr(GOFChi)",
                  "Test", "   Df", "LR stat.", "Pr(Chi)")
  if(test=="none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of Loglinear Models\n")
  out
}

