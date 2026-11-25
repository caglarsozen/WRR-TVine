# Fixed-parameter GARCH filter causality and marginal density/inverse checks.
# Synthetic data only; no empirical accuracy or convergence guarantee is implied.
source("R/wrr_tvine.R", local = TRUE, encoding = "UTF-8")
stopifnot(requireNamespace("rugarch",quietly=TRUE))
set.seed(7021)
n<-650L; e<-rt(n,df=7)/sqrt(7/5); x<-numeric(n); h<-rep(.0001,n)
for(i in seq_len(n)) {
  if(i>1) h[i]<-.000005+.08*x[i-1]^2+.87*h[i-1]
  x[i]<-sqrt(h[i])*e[i]
}
train<-seq_len(n)<=400L; valid<-seq_len(n)>400L & seq_len(n)<=500L
for(model in c("sGARCH","gjrGARCH","eGARCH")) {
  m<-wrr_with_seed(7031,wrr_garch_margin(x,train,model))
  stopifnot(all(is.finite(m$u)),all(m$sigma>0),is.finite(wrr_margin_nll(m,x,valid)))
  back<-wrr_margin_inverse(m,m$u[501:550],501:550)
  stopifnot(max(abs(back-x[501:550]))<1e-7)
  # Fixed training parameters; perturb future observations after date 550.
  changed<-x;changed[551:n]<-20*changed[551:n]+.1
  f<-rugarch::ugarchfilter(m$fixed,data=changed,n.old=sum(train),rec.init="all")
  stopifnot(max(abs(as.numeric(rugarch::fitted(f))[1:550]-m$mu[1:550]))<1e-10)
  stopifnot(max(abs(as.numeric(rugarch::sigma(f))[1:550]-m$sigma[1:550]))<1e-10)
  prefix<-rugarch::ugarchfilter(m$fixed,data=x[1:550],n.old=sum(train),rec.init="all")
  stopifnot(max(abs(as.numeric(rugarch::sigma(prefix))-m$sigma[1:550]))<1e-10)
  cat("PASS",model,"density, inverse-CDF, future-perturbation and prefix checks\n")
}
cat("GARCH INTERFACE AND CAUSALITY CHECKS PASSED.\n")
