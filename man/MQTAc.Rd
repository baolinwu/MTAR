% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlm.R
\name{MQTAc}
\alias{MQTAc}
\title{Multiple quantitative trait association test with common covariates}
\usage{
MQTAc(obj, G)
}
\arguments{
\item{obj}{fitted null model from MLM.null}

\item{G}{genotype vector}
}
\value{
\describe{
  \item{p.value}{ three association p-values: an ombinus m-DF Wald test; two 1-DF Wald tests assuming common effect or common scaled effect (see ref) }
  \item{coef}{ estimated variant regression coefficients for all traits}
}
}
\description{
Extremely efficient computation of genome-wide association test of multiple quantitative traits.
}
\references{
Wu,B. and Pankow,J.S. (2018) Fast and accurate genome-wide association test of multiple quantitative traits. \emph{Computational and mathematical methods in medicine}, in press.
}
