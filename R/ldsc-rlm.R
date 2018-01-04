#' A robust and efficient LD score regression (LDSC) for estimating the SNP heritability using GWAS summary data
#'
#' Follow the Bulik-Sullivan et al. (2015) approach but use a weighted robust LS estimation without filtering large summary statistics.
#' Here the input r2 is the averaged reference LD scores (scaled by the total number of SNPs used to compute the LD scores).
#' The weight W is typically based on the LD scores computed using HapMap3 common SNPs. The purpose is to correct over-counting.
#' LDSC is a variance regression (VR) approach: regressing chi-square statistics on LD scores.
#'
#' @param Z summary Z-statistics for M variants
#' @param r2 average reference LD scores for M variants
#' @param N GWAS sample size for each variant (could be different across variants)
#' @param W variant weight
#'
#' @export
#' @return
#' \describe{
#'   \item{h2}{ SNP heritability }
#'   \item{v0}{ intercept in the LDSC, related to the genomic control (GC) parameter }
#' }
#' @references
#' Bulik-Sullivan,B.K. et al. (2015) LD Score regression distinguishes confounding from polygenicity in genome-wide association
#' studies. Nature Genetics, 47, 291-295.
#'
#' Guo,B. and Wu,B. (2017) Principal component based adaptive association test of multiple traits using GWAS summary statistics.
SHvr <- function(Z,r2, N, W=NULL){
  if(is.null(W)) W = rep(1,length(Z))
  tau = (mean(Z^2)-1)/mean(N*r2)
  rcf = as.vector( rlm(I(Z^2) ~ I(N*r2), weight=W/(1+tau*N*r2)^2, psi=Ghuber)$coef )
  return(list(h2=rcf[2], v0=rcf[1]) )
}




#' A robust and efficient LD score regression approach to estimating the genetic correlation using GWAS summary data
#'
#' Follow the Bulik-Sullivan et al. (2015) approach but use a weighted robust LS estimation. 
#' Here the input r2 is the averaged reference LD scores (scaled by the total number of SNPs used to compute the LD scores).
#' The weight W is typically based on the LD scores computed using HapMap3 common SNPs to correct over-counting. 
#'
#' @param Z Mx2 matrix of summary Z-statistics for M variants from two GWAS
#' @param r2 average reference LD scores for M variants
#' @param N1 sample size for the 1st GWAS
#' @param N2 sample size for the 2nd GWAS
#' @param Nc overlapped sample size between two GWAS
#' @param W variant weight
#' 
#' @return
#' \describe{
#'   \item{gv}{ genetic covariance }
#'   \item{gc}{ genetic correlation }
#'   \item{r0}{ estimated intercept, related to the marginal trait correlation }
#'   \item{h2s}{ SNP heritabilities }
#' }
#'
#' @export
#' @references
#' Bulik-Sullivan,B.K. et al. (2015) An atlas of genetic correlations across human diseases and traits. Nat. Genet. 47, 1236–1241.
#'
#' Guo,B. and Wu,B. (2017) Principal component based adaptive association test of multiple traits using GWAS summary statistics. tech rep.
GCvr <- function(Z,r2, N1,N2,Nc=0, W=NULL){
  if(is.null(W)) W = rep(1,length(r2))
  h1 = SHvr(Z[,1],r2, N1, W)$h2
  h2 = SHvr(Z[,2],r2, N2, W)$h2
  Y = Z[,1]*Z[,2]
  X = (sqrt(N1)*sqrt(N2)+sqrt(Nc/N1*N2))*r2
  N1r2 = N1*r2; N2r2 = N2*r2
  ## 1st round
  if(any(Nc>0)){
    rcf = as.vector(rlm(Y ~ X, psi=Ghuber)$coef)
    intc = rcf[1]; gv = rcf[-1]
  } else{
    intc = 0
    gv = as.vector(rlm(Y ~ X-1, psi=Ghuber)$coef)
  }
  ## 2nd round
  Wt = W/( (h1*N1r2+1)*(h2*N2r2+1) + (X*gv + intc)^2 )
  if(any(Nc>0)){
    rcf = as.vector( rlm(Y ~ X, weight=Wt, psi=Ghuber)$coef )
    intc = rcf[1]; gv = rcf[-1]
  } else{
    intc = 0
    gv = as.vector(rlm(Y ~ X-1, weight=Wt, psi=Ghuber)$coef)
  }
  gc = gv/sqrt(h1*h2)
  return(list(gv=gv,gc=gc, r0=intc, h2s = c(h1,h2)) )
}


