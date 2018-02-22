KATqval = function(pval, lambda, racc=1e-2,lim=1e7){
  ## safe check
  if( all(abs(lambda-lambda[1])/max(abs(lambda))<1e-10) ){
    return( qchisq(pval,length(lambda),lower=FALSE)*mean(lambda) )
  }
  piv = floor(1/pval)
  f0 = function(x) ( KATpval(x,lambda,racc,lim) - pval )*piv
  q0 = Liu.qval(pval, lambda)
  a1 = q0; a2 = q0*1.5
  while(f0(a1)*f0(a2)>0){
    a1 = a1/2; a2 = a2*2
  }
  uniroot(f0, c(a1,a2), tol = q0*1e-4)$root
}
Liu.qval = function(pval, lambda){
  param = Liu.param(lambda)
  Qx = qchisq(pval,df=param$l,ncp=param$d, lower.tail=FALSE)
  (Qx - param$muX)/param$sigmaX*param$sigmaQ + param$muQ
}

