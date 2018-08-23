#' Various multi-trait association tests using GWAS summary data
#'
#' We compute analytical p-values for three multi-trait tests: minimum marginal test p-value (minP), sum of Z-stats (ST),
#' and sum of squared Z-stats (SQT). 
#'
#' @param  Z summary Z-statistics across multiple traits
#' @param  Sig the estimated marginal trait correlation matrix
#' @return vector of p-values for three tests:  minP,ST,SQT
#' @export
#' @references
#' Massoti M*, Guo B* and Wu B. (2018) Pleiotropy informed adaptive association test of multiple traits using GWAS summary data. \emph{Biometrics}, under revision (* contribute equallly).
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
  p.value = c(pvalm,pvals,pval2)
  names(p.value) = c('minP', 'ST', 'SQT')
  return(p.value)
}
