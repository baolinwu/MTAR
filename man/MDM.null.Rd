% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdm.R
\name{MDM.null}
\alias{MDM.null}
\title{Null multivariate linear regression model fitting}
\usage{
MDM.null(YX, pux = NULL)
}
\arguments{
\item{YX}{list of outcome and covariates: first component is matrix of outcomes (samples in rows)}

\item{pux}{the number of covariates counted when calcuating degree of freedoms}
}
\value{
\describe{
  \item{vlm}{ fitted null model }
  \item{n,m,pux}{ sample size and outcome/covariate dimensions }
  \item{f0}{ model formula}
  \item{xlist}{ constraint list}
  \item{YX}{ data list}
}
}
\description{
Outcomes are regressed on different set of covariates.
}
