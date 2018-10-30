#' Principal component based adaptive multi-trait association test using GWAS summary data
#'
#' We show that using the summary statistics, we can compute the PC to summarize the multiple traits,
#' and then construct the 1-DF test statistic (ET). The omnibus K-DF chi-square test (OT) is generally robust and powerful.
#' We then define the adaptive test (AT) as the minium p-values of weighted sums of ET (\eqn{\rho}) and OT (\eqn{1-\rho}) tests.
#' Efficient algorithms are developed to compute the analytical p-values for all tests.
#' We use the LD score regression (see GCvr() function) to accurately estimate the marginal trait correlation using GWAS summary data.
#'
#' @param  Z summary Z-statistics across multiple traits
#' @param  Sig the estimated marginal trait correlation matrix
#' @param  rho sequence of weights assigned to the 1-DF test
#' @return
#' \describe{
#'   \item{p.value}{ p-value for the AT }
#'   \item{pvals}{ the list of all p-values }
#'   \item{rho.est}{ estimated optimal \eqn{\rho} value }
#' }
#' @export
#' @references
#' Bulik-Sullivan B et al. (2015) An atlas of genetic correlations across human diseases and traits. \emph{Nature Genetics}, 47(11):1236--41.
#' 
#' Guo,B. and Wu,B. (2018) Principal component based adaptive association test of multiple traits using GWAS summary statistics. bioRxiv 269597; doi: 10.1101/269597
#' 
#' Guo,B. and Wu,B. (2018) Integrate multiple traits to detect novel disease-gene association using GWAS summary data with an adaptive test approach. \emph{Bioinformatics}, under revision.
emats <- function(Z,Sig, rho=0:5/5){
  K = length(rho); M = dim(Sig)[1]
  es = eigen(Sig,sym=TRUE)
  chi1 = sum(Z*solve(Sig,Z))
  chi2 = sum(es$vec[,1]*Z)^2/es$val[1]
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
    return(list(p.value=0, pvals=pval, rho.est=rho[which.min(pval)]) )
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
  return(list(p.value=p.val, pvals=pval, rho.est=rho[which.min(pval)]) )
}



