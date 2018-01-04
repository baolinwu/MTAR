#' Compute the tail probability of weighted sum of 1-DF chi-square rvs
#'
#' Use Davies' method and apply a relative error bound (avoid redundant computations): efficient and accurate.
#' @param Q.all  test statistics
#' @param lambda  mixing coefficients
#' @param acc  relative error bound
#' @param lim  maximum number of integration terms
#' @return tail probability
#' @export
#' @references
#' Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. \emph{Annals of human genetics}, 80(2), 123-135.
#'
#' Guo,B. and Wu,B. (2017) Principal component based adaptive association test of multiple traits using GWAS summary statistics. tech rep.
KATpval <- function(Q.all, lambda, acc=1e-2,lim=1e7){
  pval = Liu.pval(Q.all,lambda)
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = Wdavies(Q.all[i],lambda,acc=acc*pval[i],lim=lim); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}



### saddlepoint approx: modified from Lumley survey package.
saddle = function(x,lambda){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0]))*0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
  }
}
Sadd.pval = function(Q.all,lambda){
  sad = rep(1,length(Q.all))
  if(sum(Q.all>0)>0){
    sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda)
  }
  id = which(is.na(sad))
  if(length(id)>0){
    sad[id] = Liu.pval(Q.all[id], lambda)
  }
  return(sad)
}
### modified Liu method; internal function
Liu.param = function(lambda){
    c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
    muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
    s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
    if(s1^2 > s2){
      a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
    } else {
      l = 1/s2;  a = sqrt(l);  d = 0
    }
    muX = l+d;  sigmaX = sqrt(2)*a
    list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
}
LiuPval = function(Q.all, param){
  Qx = (Q.all - param$muQ)/param$sigmaQ*param$sigmaX + param$muX
  pchisq(Qx, df = param$l,ncp=param$d, lower.tail=FALSE)
}
Liu.pval = function(Q, lambda){
  param = Liu.param(lambda)
  LiuPval(Q, param)
}
Liu0.qval = function(pval, lambda){
  param = Liu.param(lambda)
  df = param$l
  Qx = qchisq(pval,df=df,lower.tail=FALSE)
  (Qx - df)/sqrt(2*df)*param$sigmaQ + param$muQ
}
