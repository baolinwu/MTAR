#' Null multivariate linear regression model fitting
#'
#' All outcomes are regressed on the same set of covariates.
#' @param  Ys  matrix of outcomes (samples in rows and trait in columns)
#' @param  Xs  matrix of common covariates for all traits. Default to NULL for no covariates.
#' @return
#' \describe{
#'   \item{res}{ residual matrix (sample in rows) }
#'   \item{Ux}{ eigen vectors of projection matrix }
#'   \item{n,m,p}{ sample size and outcome/covariate dimensions }
#' }
#' @export
## null model fitting based on SVD
MLM.null <- function(Ys,Xs=NULL){
  n = dim(Ys)[1]; m = dim(Ys)[2]
  if(is.null(Xs)){
    p = 0
    Ux = matrix(1/sqrt(n),n,1)
  } else{
    p = dim(Xs)[2]
    X = cbind(1,Xs)
    ## SVD
    Ux = svd(X,nv=0)$u
  }
  res = Ys - Ux%*%(t(Ux)%*%Ys)
  return(list(res=res,Ux=Ux,n=n,m=m,p=p))
}

#' Multiple quantitative trait association test with common covariates
#'
#' Extremely efficient computation of genome-wide association test of multiple quantitative traits.
#' @param  obj fitted null model from MLM.null
#' @param  G genotype vector
#' @return
#' \describe{
#'   \item{p.value}{ three association p-values: an ombinus m-DF chi-square test; two 1-DF chi-square tests assuming common effect or common scaled effect (see ref) }
#'   \item{coef}{ estimated variant regression coefficients for all traits}
#' }
#' 
#' @export
#' @references
#' Wu,B. and Pankow,J.S. (2017) Fast and accurate genome-wide association test of multiple quantitative traits. tech report.
MQTAc <- function(obj,G){
  n = obj$n; m=obj$m; p=obj$p
  pval = rep(NA,3)
  Gh = as.vector( obj$Ux%*%colSums(obj$Ux*G) )
  Ge = G-Gh; G2 = sum(Ge^2)
  if(G2<1e-3) return(list(p.value=pval,coef=NA))
  ## beta coef
  U0 = colSums(obj$res*Ge)
  rcf = U0/G2
  res = obj$res - outer(Ge,rcf)
  R0 = t(res)%*%res/(n-p-2); R0i = solve(R0); sigs = sqrt(diag(R0))
  ## omnibus test
  RU0 = R0i%*%U0
  Z = sum(U0*RU0)/G2
  Zf = Z*(n-p-m-1)/m/(n-p-2)
  pval[1] = pf(Zf,m,n-p-m-1,lower=FALSE)
  ## 1-DF common scaled effect
  Z2 = sum(RU0*sigs)^2/G2/sum(t(R0i*sigs)*sigs)
  Zf2 = Z2*(n-p-m-1)/(n-p-2)
  pval[2] = pf(Zf2,1,n-p-m-1,lower=FALSE)
  ## 1-DF common beta
  Z1 = sum(RU0)^2/G2/sum(R0i)
  Zf1 = Z1*(n-p-m-1)/(n-p-2)
  pval[3] = pf(Zf1,1,n-p-m-1,lower=FALSE)
  ## return(list(p.value=pval,coef=rcf,Sigma=R0))  
  return(list(p.value=pval,coef=rcf))
}

#' Minimum p-value based multiple quantitative trait association test
#'
#' We efficiently computate minimum test p-values (minP) across multiple quantitative traits;
#' and report two test p-values: 1) Bonferroni corrected minimum p-values (Pbonf) across traits.
#' 2) Analytical significance p-value (Pmin) of minP using asymptotic multivariate normal integration.
#' See refs in the PGmvn() function.
#'
#' Remarks: (1) Be cautious to interpret extreme Pmin: it is very hard to accurately compute extreme p-values.
#' (2) Generally Pmin is close to Pbonf at extreme values under moderate trait correlations.
#' (3) Under extreme trait correlations, Pmin can offer advantages compared to Pbonf.
#' (4) Note that theoretically we have minP \eqn{\le} Pmin \eqn{\le} Pbonf.
#' @param  obj fitted null model from MLM.null
#' @param  G genotype vector
#' @return
#' \describe{
#'   \item{p.value}{ individual trait association test p-values }
#'   \item{Pbonf}{ Bonferroni corrected minimum p-value }
#'   \item{Pmin}{ Significance p-value of minimum p-value computed analytically based on multivariate normal dist }
#' }
#' 
#' @export
#' @references
#' Conneely,K.N. and Boehnke,M. (2007) So many correlated tests, so little time! Rapid adjustment of P values for multiple correlated tests. Am. J. Hum. Genet. 81, 1158–1168.
#'
#' Conneely,K.N. and Boehnke,M. (2010) Meta-analysis of genetic association studies and adjustment for multiple testing of correlated SNPs and traits. Genetic Epidemiology. 34:739-746.
#' 
#' Wu,B. (2017) MTAR: an R package for versatile genome-wide association test of multiple traits. tech report.
MTA.minp <- function(obj,G){
  n = obj$n; m=obj$m; p=obj$p
  Gh = as.vector( obj$Ux%*%colSums(obj$Ux*G) )
  Ge = G-Gh; G2 = sum(Ge^2)
  if(G2<1e-3) return(list(p.value=rep(NA,m),coef=NA))
  ## beta coef
  U0 = colSums(obj$res*Ge)
  rcf = U0/G2
  res = obj$res - outer(Ge,rcf)
  S2 = colSums(res^2)/(n-p-2)
  sigs = sqrt(S2)
  Z = rcf/sigs*sqrt(G2)
  pval = pt(-abs(Z),df=n-p-2)*2
  ## uni-test
  R0 = t(res)%*%res/(n-p-2)
  R0s = t(R0/sigs)/sigs
  p0 = min(pval)
  pbonf = min(p0*m, 1)
  Z0 = -qnorm(p0/2)
  psvd = try( {PGmvn(lower=rep(-Z0,m),upper=rep(Z0,m),sigma=R0s,Nsample=1e5,Ncov=m)} )
  if(class(psvd)=='try-error'){
    it = 0; maxpts=1e5;  ## abseps = max(pbonf*0.05,1e-8);  ## abseps = min(max(pbonf*0.05, p0/2), 1e-3)
    abseps = min(max(pbonf*0.05, p0), 1e-5)
    pmvn = 1-pmvnorm(lower=rep(-Z0,m),upper=rep(Z0,m), sigma=R0s, algorithm=GenzBretz(maxpts=maxpts,abseps=abseps))
    while( (attr(pmvn,'msg')!='Normal Completion')&(it<10) ){
      it = it+1;  maxpts = maxpts*2
      pmvn = 1-pmvnorm(lower=rep(-Z0,m),upper=rep(Z0,m), sigma=R0s, algorithm=GenzBretz(maxpts=maxpts,abseps=abseps))
    }
    pmin = min(max(pmvn, p0), pbonf)
  } else{
    pmin = min(max(psvd, p0), pbonf)
  }
  return(list(p.value=pval, pmin=pmin,pbonf=pbonf))
  ## Genz() is not accurate!!! (quasi-MC) almost not working for extreme p-values!!!
  ## Miwa() can produce negative values!!! and computationally intensive for large m!!!
  ## PGmvn() can often lead to much more accurate results.
  
  ## return(list(p.value=pval, pmin=pmin,pbonf=pbonf, psvd=psvd))
  ## it = 0; steps=128
  ## pmiwa = 1-pmvnorm(lower=rep(-Z0,m),upper=rep(Z0,m), sigma=R0s, algorithm=Miwa(steps=steps))
  ## while( (attr(pmiwa,'msg')!='Normal Completion')&(it<10) ){
  ##  it = it+1;  steps = steps*2
  ##  pmiwa = 1-pmvnorm(lower=rep(-Z0,m),upper=rep(Z0,m), sigma=R0s, algorithm=Miwa(steps=steps))
  ## }
  ## return(list(p.value=pval, pmin=pmin,pbonf=pbonf, pmiwa=pmiwa))
  ## 1-psvdnorm(lower=rep(-Z0,m),upper=rep(Z0,m),Sigma=R0s,Nsample=1e5,Ncov=m)
}




#' Compute the tail probability of the m-dim multivariate normal distribution
#'
#' Internal function. Not to be called directly.
#' @param  lower the vector of lower limits of length m
#' @param  upper the vector of upper limits of length m
#' @param  mean the mean vector of length m
#' @param  sigma the covariance matrix of dimension m
#' @param  Nsample the number of Monte Carlo samples
#' @param  Ncov the number of control variates to be used (<=m).
#' @return multivariate normal distribution probability of outside the specified box region.
##' @export
#' @references
#' Phinikettos,I. and Gandy,A. (2011) Fast computation of high-dimensional multivariate normal probabilities. Computational Statistics & Data Analysis. 55, 1521–1529.
#'
#' Genz, A., Bretz, F., Miwa, T., Mi, X., Leisch, F., Scheipl, F., Bornkamp, B., Maechler, M., Hothorn, T. (2015) mvtnorm: Multivariate Normal and t Distributions. R package version 1.0-3. \url{https://cran.r-project.org/web/packages/mvtnorm/index.html}
PGmvn <- function(lower=-Inf, upper=Inf, mean=NULL, sigma, Nsample=1e4,Ncov=1){
  m = dim(sigma)[1]
  if(is.null(mean)){
    mu = rep(0,m)
  } else{
    mu = mean
  }
  if(length(lower)<=1) lower = rep(lower,m)
  if(length(upper)<=1) upper = rep(upper,m)
  lz = lower-mu; uz = upper-mu
  1-psvdnorm(lz,uz,Sigma=sigma,Nsample=Nsample,Ncov=Ncov)
}

