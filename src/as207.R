#
# Iterative proportional fitting EM approach to fitting log-linear
# models to incomplete tables.
# Modelled after Algorithm AS 207 Appl. Statist. (1984) vol.33, no.3
#
emgllm <- function(y,s,X,maxit=1000,tol=0.00001) {
  eps<-1e-5
  if (maxit<=0) maxit<-1
  if (typeof(X)=="language") {
    X<-model.matrix(X)
  }
  X<-cbind(X,double(nrow(X)))
  nkk<-nrow(X)
  full.table<-rep(1.0,nrow(X))
#
# standardize the design matrix
  X<-X-min(X)
  while(abs((maxrs<-max(rs<-rowsum(X))-1.0)>eps) {
    X<-X/maxrs
  }
  if (any(abs(rs-1.0)>eps) {
    X[,nkk]<-1.0-rs
  }
#
# enter the EM algorithm
  it<-0
  while(TRUE) {
    it<-it+1
    f<-s %*% full.table
    if (it>itmax || all(abs((f-oldf)<tol))) {
      break
    }
    oldf<-f
    w<-y[s]
    idx<-f>eps
    w[idx]<-full.table[idx] %*% y[idx] / f[idx]
    v1 <- X*w
#
# enter the IPF algorithm
    while (all(abs(full.table-w2)>conv)) {
      w2 <- full.table
     
      v2 <- double(ncol(X))
      v2 <- v2 + X(i,k) * full.table(i)
      idx<- v2>eps
     
      full.table[idx] <- (full.table[idx] * v1/v2) ^ X(i, k)
    }
  }
  idx<-(y>0 & f>0)
  deviance <- 2*sum(y[idx] * log(y[idx] / f[idx]))
  list(deviance=deviance,
       observed.values=y,
       fitted.values=f,
       full.table=full.table)
}
