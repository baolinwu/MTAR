## Function to compute AT p-value as a function of the minimum p-value
Tval = function(minp, Sig, rho=0:5/5){
  K = length(rho); M = dim(Sig)[1]
  qval = rep(0,K)
  for(i in 1:K){
    if(rho[i]==0){
      qval[i] = qchisq(minp,M, lower=FALSE)
    } else if(rho[i]==1){
      qval[i] = qchisq(minp,1, lower=FALSE)
    } else{
      qval[i] = KATqval(minp, c(rep(1-rho[i],M-1),1))
    }
  }
  chiint = function(xq){
    xi = sapply(xq, function(x) min( (qval[-K]-x)/(1-rho[-K])) )
    pchisq(xi,M-1,lower=FALSE)*dchisq(xq,1)
  }
  p.val = minp + integrate(chiint, 0,qval[K])$val
  attr(p.val,'qval')  = qval
  return(p.val)
}


#' Compute the list of significant SNPs using the GWAS summary data based multi-trait association tests
#'
#' We compute the significant SNPs for three tests: omnibus chi-square test (OT); 1-DF PC test (ET); adaptive test (AT).
#'
#' @param  Z matrix of summary Z-statistics (SNPs by traits)
#' @param  Sig the estimated marginal trait correlation matrix
#' @param  rho sequence of weights assigned to the 1-DF test
#' @param  alpha desired genome-wide significance level (default to 5E-8)
#' @return
#' \describe{
#'   \item{idQ}{ significant SNP list for OT }
#'   \item{idE}{ significant SNP list for ET }
#'   \item{idA}{ significant SNP list for AT }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2018) Integrate multiple traits to detect novel disease-gene association using GWAS summary data with an adaptive test approach. \emph{Bioinformatics}, under revision.
Lemats <- function(Z, Sig, rho=0:5/5, alpha=5e-8){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Sig)[1]; K = length(rho)
  t1 = qchisq(alpha,M, lower=FALSE)
  t2 = qchisq(alpha,1, lower=FALSE)
  ##
  f2 = function(minp) (Tval(minp,Sig,rho) - alpha)*floor(1/alpha)
  p0 = uniroot(f2, c(alpha/K,alpha), tol=alpha*1e-4)$root
  aa = Tval(p0,Sig,rho)
  Qs = attributes(aa)$qval
  ##
  Sigi = solve(Sig)
  es = eigen(Sig, sym=TRUE)
  Q = rowSums(Z%*%Sigi*Z)
  B = (Z%*%es$vec[,1])^2/es$val[1]
  ss1 = which(Q>t1)
  ss2 = which(B>t2)
  ss3 = 1*( ((1-rho[1])*Q+rho[1]*B)>Qs[1] )
  for(i in 2:K)  ss3 = ss3 + 1*( ((1-rho[i])*Q+rho[i]*B)>Qs[i] )

  return(list(idQ=ss1,idE=ss2,idA=which(ss3>0)) )
}


#' Compute the list of significant SNPs using the GWAS summary data based multi-trait association tests
#'
#' We compute the significant SNPs for (1) minimum marginal test p-value (minP); SZ test; SZ2 test;
#' and (2) omnibus chi-square test (OT); minimum of OT and SZ2 (MU); and significant at any marginal test.
#'
#' @param  Z matrix of summary Z-statistics (SNPs by traits)
#' @param  Sig the estimated marginal trait correlation matrix
#' @param  alpha desired genome-wide significance level (default to 5E-8)
#' @return
#' \describe{
#'   \item{idM}{ significant SNP list for minP }
#'   \item{idSZ}{ significant SNP list for SZ }
#'   \item{idSZ2}{ significant SNP list for SZ2 }
#'   \item{idOT}{ significant SNP list for OT }
#'   \item{idMU}{ significant SNP list for min(SZ2,OT) }
#'   \item{idm}{ SNPs significant at any marginal trait }
#' }
#' @export
#' @references
#' Guo,B. and Wu,B. (2018) Integrate multiple traits to detect novel disease-gene association using GWAS summary data with an adaptive test approach. \emph{Bioinformatics}, under revision.
Lmatz <- function(Z, Sig, alpha=5e-8){
  if(class(Z)!='matrix') Z = as.matrix(Z)
  M = dim(Sig)[1]
  es = eigen(Sig, sym=TRUE)
  ##
  ## f1 = function(q0) ( 1-pmvnorm(rep(-q0,M),rep(q0,M), sigma=Sig, algorithm=Miwa(steps=256)) - alpha )*floor(1/alpha)
  f1 = function(q0) ( 1-pmvnorm(rep(-q0,M),rep(q0,M), sigma=Sig) - alpha )*floor(1/alpha)
  q1 = qnorm(alpha/2,lower=FALSE); q2 = q1*1.1
  while( f1(q1)*f1(q2)>0 ) q2 = q2*1.1
  q0 = uniroot(f1, c(q1,q2), tol=alpha*1e-4)$root
  ##
  t1 = qchisq(alpha,1, lower=FALSE)
  t2 = KATqval(alpha,es$val)
  ##
  SZ1 = rowSums(Z)^2/sum(Sig)
  ss1 = which(SZ1>t1)
  SZ2 = rowSums(Z^2)
  ss2 = which(SZ2>t2)
  ## Zm = abs(Z[,1]); for(i in 2:M)  Zm = pmax(Zm, abs(Z[,i])); ss3 = which(Zm>q0)
  ss3 = which(rowSums(abs(Z)>q0)>0)
  ss0 = which(rowSums(abs(Z)>q1)>0)
  ##
  tm = qchisq(alpha,M, lower=FALSE)
  Q = rowSums(Z%*%solve(Sig)*Z)
  ss4 = which(Q>tm)
  ##
  ## tm2 = qchisq(alpha/2,M, lower=FALSE)
  ## t22 = KATqval(alpha/2,es$val)
  ## ss5 = which( (Q>tm2)|(SZ2>t22) )
  t1 = qchisq(alpha/3,M, lower=FALSE)
  t0 = KATqval(alpha/3,es$val)
  th = KATqval(alpha/3,es$val/2+0.5)
  ss5 = which( (Q>t1)|(SZ2>t0)|((Q/2+SZ2/2)>th) )
  ##
  return(list(idM=ss3,idSZ=ss1,idSZ2=ss2, idOT=ss4,idMU=ss5, idm=ss0) )
}


