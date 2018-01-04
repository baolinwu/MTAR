## internal func.
# ' Compute the tail prob for weighted sum of 1-DF chi-square random variables
# '
# ' Imported from CompQuadForm R package (version 1.4.3; 10/4/2017). Slight simplifications for sum of central chi-squares.
# ' @param qx quantile to compute tail prob
# ' @param lambda weights for 1-DF chi-square rvs
# ' @param delta non-centrality parameters for 1-DF chi-square rvs. Default to zeros.
# ' @param lim total number of integration terms
# ' @param acc accuracy for computed probability
# ' 
# ' @return
# ' \describe{
# '  \item{Qq}{ tail prob }
# '  \item{trace}{vector, indicating performance of procedure, with the following components: 1: absolute value sum, 2: total number of integration terms, 3: number of integrations, 4: integration interval in main integration, 5: truncation point in initial integration, 6: standard deviation of convergence factor term, 7: number of cycles to locate integration parameters}
# '  \item{ifault}{ fault indicator: 0: no error, 1: requested accuracy could not be obtained, 2: round-off error possibly significant, 3: invalid parameters, 4: unable to locate integration parameters }
# ' }
# ' @export
Wdavies <- function(qx, lambda, delta=NULL, lim=1e9, acc=1e-6){
  m = length(lambda)
  if(is.null(delta)) delta = rep(0,m)
  h = rep(1,m); sigma = 0
  out = .C("qfc", as.double(lambda), as.double(delta),  as.integer(h),
           as.integer(m), as.double(sigma), as.double(qx), as.integer(lim),
           as.double(acc), trace = as.double(rep(0, 7)), ifault = as.integer(0),
           res = as.double(0), PACKAGE = "MTAR")
  out$res = 1 - out$res
  return(list(Qq=out$res, trace=out$trace, ifault=out$ifault))
  
}


#' @useDynLib MTAR
NULL
#> NULL

.onUnload <- function (libpath) {
  library.dynam.unload("MTAR", libpath)
}



