\name{makeBootCIs}
\alias{makeBootCIs}
\title{Calculate percentile confidence intervals from a matrix of bootstrapped predictions}
\usage{
  makeBootCIs(preds, quants = c(0.025, 0.975))
}
\arguments{
  \item{preds}{matrix of bootstrap predictions where each
  column is a bootstrap realisation}

  \item{quants}{(\code{default = c(0.025, 0.975)}. Vector
  of length two of quantiles.}
}
\description{
  Calculate percentile confidence intervals from a matrix
  of bootstrapped predictions
}
\examples{
\dontrun{
makeBootCIs(bootPreds)
}
}

