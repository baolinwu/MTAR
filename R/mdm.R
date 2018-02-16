#' Null multivariate linear regression model fitting
#'
#' Outcomes are regressed on different set of covariates.
#' @param  YX  list of outcome and covariates: first component is matrix of outcomes (samples in rows)
#' @param  pux  the number of covariates counted when calcuating degree of freedoms
#' @return
#' \describe{
#'   \item{vlm}{ fitted null model }
#'   \item{n,m,pux}{ sample size and outcome/covariate dimensions }
#'   \item{f0}{ model formula}
#'   \item{xlist}{ constraint list}
#'   \item{YX}{ data list}
#' }
#' 
#' @export
## null model fitting based on SVD
MDM.null = function(YX,pux=NULL){
### YX: a list of Ys and Xs
  names(YX)[1] = 'SurY'
  n = dim(YX$SurY)[1]; m = dim(YX$SurY)[2]
  for(k in 1:m) names(YX)[1+1:m] = paste0('SurX',1:m)
  ## intercept only counted once
  if(is.null(pux)) pux = max(sapply(YX[-1], dim)[2,])+1
  cm1 = diag(m)
  xlist = vector('list',m+1)
  xlist[[1]] = cm1; names(xlist)[[1]] = '(Intercept)'
  xlist[[2]] = cm1[,1,drop=FALSE]
  names(xlist)[2] = 'SurX1'
  f0 = paste0('SurY~SurX1')
  for(k in 2:m){
    f0 = paste0(f0, '+SurX', k)
    xlist[[k+1]] = cm1[,k,drop=FALSE]
    names(xlist)[k+1] = paste0('SurX',k)
  }
  ## SUR mle
  vlm0 = vglm(as.formula(f0), SURff(mle.normal=TRUE), data=YX, constraints=xlist, epsilon=1e-12, maxit=1e2, stepsize=0.5)
  return(list(vlm=vlm0,n=n,m=m,pux=pux,f0=f0, xlist=xlist,YX=YX))
}

#' Multiple quantitative trait association test with differing covariates
#' 
#' Fast computation of genome-wide association test of multiple quantitative traits.
#' @param  obj fitted null model from MDM.null
#' @param  G genotype vector
#' @return
#' \describe{
#'   \item{p.value}{ three association p-values: an ombinus m-DF Wald test; two 1-DF Wald tests assuming common effect or common scaled effect (see ref) }
#'   \item{coef}{ estimated variant regression coefficients for all traits}
#' }
#' 
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2018) Fast and accurate genome-wide association test of multiple quantitative traits. \emph{Computational and mathematical methods in medicine}, in press.
#' @export
MQTAd = function(obj,G){
  n = obj$n; m = obj$m; pux = obj$pux
  xlist = obj$xlist; xlist[[m+2]] = diag(m); names(xlist)[m+2] = 'G'
  f1 = paste0(obj$f0, '+G')
  ## SUR mle
  rcf0 = coef(obj$vlm)
  vlm1 = vglm(as.formula(f1), SURff(mle.normal=TRUE), data=obj$YX, constraints=xlist, coefstart=c(rcf0,rep(0,m)), epsilon=1e-12, maxit=1e2)
  ## Wald stat
  jb1 = length(rcf0)+1:m
  rcf1 = coef(vlm1)[jb1]
  V = vcov(vlm1)[jb1,jb1];  Vi = solve(V)
  pval = rep(NA,3)
  ## omnibus test
  SDF = n-pux-m
  U = Vi%*%rcf1
  Q = sum(rcf1*U)
  pval[1] = pf(Q/m*SDF/n, m,SDF,lower=FALSE)
  ## 1-DF common beta
  Z1 = sum(U)^2/sum(Vi)
  Zf1 = Z1*SDF/n
  pval[3] = pf(Zf1,1,SDF,lower=FALSE)
  ## 1-DF common scaled effect
  sigs = sqrt(diag(V))
  Z2 = sum(U*sigs)^2/sum(t(Vi*sigs)*sigs)
  Zf2 = Z2*SDF/n
  pval[2] = pf(Zf2,1,SDF,lower=FALSE)
  ## return(list(p.value=pval,coef=rcf1,V=V, vlm=vlm1))
  return(list(p.value=pval,coef=rcf1))
}
