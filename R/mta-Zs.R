#' Various multi-trait association tests using GWAS summary data
#'
#' We compute analytical p-values for three multi-trait tests: minimum marginal test p-value (minP), sum of Z-stats (SZ), and sum of squared Z-stats (SZ2). 
#'
#' @param  Z summary Z-statistics across multiple traits
#' @param  Sig the estimated marginal trait correlation matrix
#' @return vector of p-values for three tests:  minP,SZ,SZ2
#' @export
#' @references
#' Guo,B. and Wu,B. (2017) Principal component based adaptive association test of multiple traits using GWAS summary statistics. tech rep.
matz <- function(Z,Sig){
  M = dim(Sig)[1]
  lam = eigen(Sig,sym=TRUE,only.val=TRUE)$val
  ## SZ
  SZ = sum(Z)^2/sum(Sig)
  pvals = pchisq(SZ,1,lower=FALSE)
  ## SZ2
  SZ2 = sum(Z^2)
  pval2 = KATpval(SZ2, lam)
  ## minp
  q0 = max(abs(Z))
  pvalm = 1-pmvnorm(rep(-q0,M),rep(q0,M), sigma=Sig, algorithm=Miwa(steps=256))
  return( c(pvalm,pvals,pval2) )
}
