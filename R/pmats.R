#' Pleiotropy informed adaptive association test across multiple traits using summary statistics
#'
#' Multi-trait association test based on GWAS summary statistics under common genetic effect assumption
#' 
#' Under pleiotropy of common effects across traits, an 1-DF pleiotropy test (PT) can be constructed.
#' The omnibus K-DF chi-square test (OT) is generally robust and powerful.
#' We then define the adaptive test (AT) as the minium p-values of weighted sums of PT (\eqn{\rho}) and OT (\eqn{1-\rho}).
#' An efficient algorithm is developed to exactly compute the analytical p-value of AT.
#' We recommend using the LD score regression (Bulik-Sullivan et al.) to accurately estimate the trait correlation matrix.
#' 
#' @param  Z summary Z-statistics across multiple traits
#' @param  Sig estimated trait correlation matrix
#' @param  rho sequence of weights assigned to the PT. Default to 0:10/10
#' @return
#' \describe{
#'   \item{p.value}{ vector of p-values for: AT, OT and PT }
#'   \item{pvals}{ the list of all p-values }
#'   \item{rho.est}{ estimated optimal \eqn{\rho} value }
#' }
#' @export
#' @references
#' Bulik-Sullivan BK  et al. (2015) LD Score regression distinguishes confounding from polygenicity in genome-wide association studies.
#'  \emph{Nature Genetics}, 47(3):291--295.
#' 
#' Bulik-Sullivan B et al. (2015) An atlas of genetic correlations across human diseases and traits. \emph{Nature Genetics}, 47(11):1236--1241.
#' 
#' Massoti M*, Guo B* and Wu B. (2018) Pleiotropy informed adaptive association test of multiple traits using GWAS summary data. \emph{Biometrics}, under revision (* contribute equallly).
#' @examples
#' Sig = diag(3)
#' Sig[c(2:3,6, 4,7,8)] = c(0.2646,-0.2317, 0.7451)
#' Z = c(1,2,4)
#' pmats(Z,Sig)
pmats <- function(Z, Sig, rho=0:10/10){
  iSig = solve(Sig)
  K = length(rho); M = dim(iSig)[1]
  tmp = iSig%*%Z;  chi1 = sum(Z*tmp);  chi2 = sum(tmp)^2/sum(iSig)
  x = (1-rho)*chi1+rho*chi2
  pval = rep(1,K)
  for(i in 1:K){
    if(rho[i]==0){
      pval[i] = pchisq(x[i],M, lower=FALSE)
    } else if(rho[i]==1){
      pval[i] = pchisq(x[i],1, lower=FALSE)
    } else{
      pval[i] = KATpval(x[i], c(rep(1-rho[i],M-1),1))
    }
  }
  minp = min(pval)
  if(minp<=0){
    p.value = c(minp, pval[rho==0], pval[rho==1])
    names(p.value) = c('AT','OT', 'PT')
    return( list(p.value=p.value, pvals=pval, rho.est=rho[which.min(pval)]) )
  }
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
  p.value = c(p.val, pval[rho==0], pval[rho==1])
  names(p.value) = c('AT','OT', 'PT')
  return(list(p.value=p.value, pvals=pval, rho.est=rho[which.min(pval)]) )
}
