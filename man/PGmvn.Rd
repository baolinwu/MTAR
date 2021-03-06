% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlm.R
\name{PGmvn}
\alias{PGmvn}
\title{Compute the tail probability of the m-dim multivariate normal distribution}
\usage{
PGmvn(lower = -Inf, upper = Inf, mean = NULL, sigma,
  Nsample = 10000, Ncov = 1)
}
\arguments{
\item{lower}{the vector of lower limits of length m}

\item{upper}{the vector of upper limits of length m}

\item{mean}{the mean vector of length m}

\item{sigma}{the covariance matrix of dimension m}

\item{Nsample}{the number of Monte Carlo samples}

\item{Ncov}{the number of control variates to be used (<=m).}
}
\value{
multivariate normal distribution probability of outside the specified box region.
}
\description{
Internal function. Not to be called directly.
}
\references{
Phinikettos,I. and Gandy,A. (2011) Fast computation of high-dimensional multivariate normal probabilities. Computational Statistics & Data Analysis. 55, 1521–1529.

Genz, A., Bretz, F., Miwa, T., Mi, X., Leisch, F., Scheipl, F., Bornkamp, B., Maechler, M., Hothorn, T. (2015) mvtnorm: Multivariate Normal and t Distributions. R package version 1.0-3. \url{https://cran.r-project.org/web/packages/mvtnorm/index.html}
}
