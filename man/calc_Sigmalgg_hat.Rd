% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em.R
\name{calc_Sigmalgg_hat}
\alias{calc_Sigmalgg_hat}
\title{Calculate Sigmalgg_hat}
\usage{
calc_Sigmalgg_hat(Xmat, Vg, Ve, Dmat, l)
}
\arguments{
\item{Xmat}{a genotypes matrix for a single marker, dimension a by n, where a is the number of allelotypes and n the number of subjects}

\item{Vg}{current value of genetic covariance matrix, with dimension d by d}

\item{Ve}{current value of error covariance matrix, with dimension d by d}

\item{Dmat}{a diagonal matrix containing the eigenvalues of the kinship matrix, K = UtDU}

\item{l}{a number from 1 to n, corresponding to the subject number}
}
\description{
Calculate Sigmalgg_hat
}
