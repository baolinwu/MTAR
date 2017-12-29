#' A multinomial regression model for testing multi-trait association
#'
#' An ACL (adjacent-category-logit) model that accounts for population stratification and tests multi-trait association.
#' @param  Ys matrix of multiple continuous traits to be tested
#' @param  Xs covariates to be adjusted
#' @param  G genotype vector
#' @return
#' \describe{
#'   \item{p.value}{ three p-values: an ombinus m-DF chi-square test; two 1-DF chi-square tests assuming common effect or common scaled effect (see ref) }
#' }
#' 
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2015) Statistical methods for association tests of multiple continuous traits in genome-wide association studies. \emph{Annals of human genetics}, 79(4), 282-293.
MTA.ACL = function(Ys,Xs,G){
  n = dim(Ys)[1]; m = dim(Ys)[2]
  pval = rep(NA,3)
  ## Qg: m-DF omnibus test
  obj0 = vglm(G~Xs, acat(parallel=TRUE))
  obj1 = vglm(G~Xs+Ys, acat(parallel=TRUE))
  pval[1] = pchisq(2*logLik(obj1)-2*logLik(obj0),m, lower=FALSE)
  ## null model fitting
  res = matrix(NA, n,m)
  for(k in 1:m){
    res[,k] = lm(Ys[,k]~Xs)$res
  }
  ## Tg: common effect
  cvm = cov(res)
  Yt2 = res%*%solve(cvm)
  Yt = rowSums(Yt2)
  obj2 = vglm(G~Xs+Yt, acat(parallel=TRUE))
  pval[2] = pchisq(2*logLik(obj2)-2*logLik(obj0),1, lower=FALSE)
  ## Tg': common scaled effect
  sigs = sqrt(diag(cvm))
  Yt = colSums(t(Yt2)*sigs)
  obj3 = vglm(G~Xs+Yt, acat(parallel=TRUE))
  pval[3] = pchisq(2*logLik(obj3)-2*logLik(obj0),1, lower=FALSE)
  return(list(p.value=pval))
}