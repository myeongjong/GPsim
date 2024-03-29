% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance.R
\name{cov_expo_iso}
\alias{cov_expo_iso}
\title{Isotropic exponential covariance function}
\usage{
cov_expo_iso(locs, covparms)
}
\arguments{
\item{locs}{A matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point in R^p}

\item{covparms}{A vector with covariance parameters in the form (variance, range)}
}
\value{
A matrix with \code{n} rows and \code{n} columns, with the (i, j) entry containing the isotropic exponenital covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
}
\description{
From a location matrix (\code{locs}) and a vector with covariance parameters (\code{covparms}), this function returns an isotropic exponential covariance matrix which is one of the simplest covariance matrices.
}
\section{Parametrization}{
 The covariance parameter vector is (variance, range) = \eqn{(\sigma^2 , r)}. The form of the covariance is \deqn{ C(x, y) = \sigma^2 exp( || x - y || / r)} where \eqn{x} and \eqn{y} are locations in R^p.
}

\examples{
# grid locations
cov_expo_iso(locs = expand.grid(c(0.25, 0.75), c(0.25, 0.75)), covparms = c(1, 0.1))

# randomly selected locations
cov_expo_iso(locs = matrix(runif(8), 4, 2), covparms = c(1, 0.1))
}
