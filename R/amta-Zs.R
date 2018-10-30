#' Adaptive multi-trait association test using GWAS summary data
#'
#' We consider a Z-stat based adaptive test (AZ) by weighting sum of Z-stats (SZ) and the sum of squared Z-stats (SZ2)
#'
#' @param  Z summary Z-statistics across multiple traits
#' @param  Sig the estimated marginal trait correlation matrix
#' @param  rho weights for the SZ
#' @return 
#' \describe{
#'   \item{p.value}{ p-values for AZ, SZ2, and SZ }
#'   \item{pval}{ the list of all p-values }
#'   \item{rho.est}{ estimated optimal \eqn{\rho} value }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2018) Integrate multiple traits to detect novel disease-gene association using GWAS summary data with an adaptive test approach. \emph{Bioinformatics}, under revision.
amatz <- function(Z,Sig, rho=0:10/10){
  M = dim(Sig)[1]
  ## Zw = Z; Rw = Sig
  eR = eigen(Sig,sym=TRUE);  lamR = abs(eR$val); eta = colSums(eR$vec)*sqrt(lamR)
  R1 = sum(eta^2); R2 = sum(eta^2*lamR)
  c2 = outer(eta,eta)
  Lamq = eigen(diag(lamR) - R2/R1^2*c2, symmetric=TRUE,only.values=TRUE)$val
  ## ST
  Qb = sum(Z)^2
  pvalb = pchisq(Qb/R1, 1,lower=FALSE)
  ## S2T
  Qv = sum(Z^2)
  pvalv = KATpval(Qv,lamR)
  ## AT
  L = length(rho)
  if(L<=2){
    return(list(p.value=c(A=NULL, S2=pvalv, S=pvalb), pval=pval) )
  } 
  L1 = L-1; rho1 = rho[-L]
  Qw = (1-rho)*Qv + rho*Qb
  pval = rep(1, L)
  pval[1] = pvalv; pval[L] = pvalb
  Lamk = vector('list', L)
  Lamk[[L]] = R1;  Lamk[[1]] = lamR
  for(k in 2:L1){
    mk = rho[k]*c2;  diag(mk) = diag(mk) + (1-rho[k])*lamR
    aak = zapsmall( abs(eigen(mk, sym=TRUE, only.val=TRUE)$val) )
    Lamk[[k]] = aak[aak>0]
    pval[k] = KATpval(Qw[k],Lamk[[k]])
  }
  minP = min(pval)
  if(minP<=0){
    return(list(p.value=c(A=0, S2=pvalv, S=pvalb), pvals=pval, rho.est=rho[which.min(pval)]) )
  }
  ##
  qval = rep(0,L1)
  for (k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  q1 = qchisq(minP,1,lower=FALSE)
  tauk = (1-rho1)*R2/R1 + rho1*R1
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    KATpval(eta1,Lamq)*dchisq(xpar,1)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
  }
  p.value = min(p.value,minP*L)
  return(list(p.value=c(A=p.value, S2=pvalv, S=pvalb), pval=pval, rho.est=rho[which.min(pval)]) )
}


## Function to compute AZ p-value as a function of the minimum p-value
Zval = function(minP, Sig, rho=0:10/10){
  eR = eigen(Sig,sym=TRUE);  lamR = abs(eR$val); eta = colSums(eR$vec)*sqrt(lamR)
  R1 = sum(eta^2); R2 = sum(eta^2*lamR)
  c2 = outer(eta,eta)
  Lamq = eigen(diag(lamR) - R2/R1^2*c2, symmetric=TRUE,only.values=TRUE)$val
  L = length(rho)
  L1 = L-1; rho1 = rho[-L]
  Lamk = vector('list', L)
  Lamk[[L]] = R1;  Lamk[[1]] = lamR
  for(k in 2:L1){
    mk = rho[k]*c2;  diag(mk) = diag(mk) + (1-rho[k])*lamR
    aak = zapsmall( abs(eigen(mk, sym=TRUE, only.val=TRUE)$val) )
    Lamk[[k]] = aak[aak>0]
  }
####
  qval = rep(0,L1)
  for (k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  q1 = qchisq(minP,1,lower=FALSE)
  tauk = (1-rho1)*R2/R1 + rho1*R1
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    KATpval(eta1,Lamq)*dchisq(xpar,1)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
  }
  p.val = min(p.value,minP*L)
  qvalL = qchisq(minP,1,lower=FALSE)*R1
  attr(p.val,'qval')  = c(qval, qvalL)
  return(p.val)
}


#' Compute the list of significant SNPs using the GWAS summary data based adaptive multi-trait association tests
#'
#' We compute the significant SNPs for the adaptive test (AZ) based on weighting the SZ and SZ2 tests.
#'
#' @param  Z matrix of summary Z-statistics (SNPs by traits)
#' @param  Sig the estimated marginal trait correlation matrix
#' @param  rho weights for the SZ
#' @param  alpha desired genome-wide significance level (default to 5E-8)
#' @return
#' \describe{
#'   \item{idAZ}{ significant SNP list for AZ }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2018) Integrate multiple traits to detect novel disease-gene association using GWAS summary data with an adaptive test approach. \emph{Bioinformatics}, under revision.
Lamatz <- function(Z, Sig, rho=0:10/10, alpha=5e-8){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  K = length(rho)
  ##
  f2 = function(minp) (Zval(minp,Sig,rho) - alpha)*floor(1/alpha)
  p0 = uniroot(f2, c(alpha/K,alpha), tol=alpha*1e-4)$root
  aa = Zval(p0,Sig,rho)
  Qs = attributes(aa)$qval
  ##
  Q = rowSums(Z^2); B = rowSums(Z)^2
  ss3 = 1*( ((1-rho[1])*Q+rho[1]*B)>Qs[1] )
  for(i in 2:K)  ss3 = ss3 + 1*( ((1-rho[i])*Q+rho[i]*B)>Qs[i] )
  return( list(idAZ=which(ss3>0)) )
}

