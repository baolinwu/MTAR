#' GEE association test of multiple quantitative traits
#'
#' An efficient implementation of genome-wide association test of multiple quantitative traits using GEE score tests.
#' @param  obj fitted null model from MLM.null
#' @param  G genotype vector
#' @return
#' \describe{
#'   \item{p.value}{ three association p-values: an ombinus m-DF chi-square test; two 1-DF chi-square tests assuming common effect or common scaled effect (see ref) }
#'   \item{U,V}{ score vector (U) and its asymptotic covariance matrix (V) }
#' }
#' 
#' @export
#' @references
#' He,Q., Avery,C.L. and Lin,D.Y. (2013) A general framework for association tests with multivariate traits in large-scale genomics studies. \emph{Genetic Epidemiology}, 37 (8), 759–767.
#'
#' Wu,B. and Pankow,J.S. (2017) Fast and accurate genome-wide association test of multiple quantitative traits. tech report.
## Marginal GEE from LinDY
MTA.lin <- function(obj,G){
  n = obj$n; m=obj$m; p=obj$p
  res = obj$res
  pval = rep(NA,3)
  ## Score test
  sigs = sqrt(colMeans(res^2))
  U = colSums(res*G)/sigs^2
  Gh = as.vector( obj$Ux%*%colSums(obj$Ux*G) )
  Ge = G-Gh; G2 = sum(Ge^2)
  if(G2<1e-3) return(list(p.value=pval,U=U,V=NA))
  Vm = t(res*Ge)/sigs^2
  V = Vm%*%t(Vm)
  ## Q: omnibus test
  Qs = sum(U*solve(V,U))
  pval[1] = pchisq(Qs,m,lower=FALSE)
  ## T: common scaled beta
  Vsig = sqrt(diag(V));  Us = U/Vsig
  R = t(V/Vsig)/Vsig; Ri = solve(R)
  Ts = sum(Ri%*%Us)/sqrt(sum(Ri))
  pval[2] = 2*pnorm(-abs(Ts))
  ## T': common beta
  Vv = diag(V); Uz = U/Vv
  R = t(V/Vv)/Vv; Ri = solve(R)
  Ts = sum(Ri%*%Uz)/sqrt(sum(Ri))
  pval[3] = 2*pnorm(-abs(Ts))
  return(list(p.value=pval,U=U,V=V))
}
