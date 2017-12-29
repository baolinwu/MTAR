# MTAR
  - an R package implementing various methods for genome-wide association test of multiple traits

-----
## Multinomial ACL model for multi-trait association test
 - References
    - Wu,B. and Pankow,J.S. (2015) Statistical methods for association tests of multiple continuous traits in genome-wide association studies. *Annals of human genetics*, 79(4), 282-293.
 - We are testing the joint effects of mm continuous traits with a SNP.
 - The current implementation requires the vglm function in the VGAM R package.
    - http://cran.r-project.org/web/packages/VGAM/index.html
 - Three tests are implemented: an omnibus m-DF chi-square test; two 1-DF chi-square tests assuming common effect or common scaled effect.


# R codes
```r
library(VGAM)
MTA.ACL = function(Ys,Xs,G){
### Ys: m continuous traits to be tested
### Xs: covariates that need to be adjusted
### G: genotype scores
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
```