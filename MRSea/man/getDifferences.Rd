\name{getDifferences}
\alias{getDifferences}
\title{Identify any significant differences between predicted data before an impact event and predicted data after an impact event}
\usage{
  getDifferences(beforePreds, afterPreds,
    quants = c(0.025, 0.975))
}
\arguments{
  \item{beforePreds}{Matrix of bootstrap predictions (n x
  B) to each grid cell before impact (same length and order
  as \code{afterPreds})}

  \item{afterPreds}{Matrix of bootstrap predictions (n x B)
  to each grid cell after impact (same length and order as
  \code{beforePreds})}

  \item{quants}{(\code{default = =c(.025,.975)}) Quantile
  for significance.}
}
\value{
  A list is returned consisting of \item{mediandiff}{Vector
  of the median difference for each grid cell}
  \item{lowerci}{Vector of the lower 2.5\% difference for
  each grid cell} \item{upperci}{Vector of the upper 97.5\%
  difference for each grid cell}
  \item{significanceMarker}{Vector of significance.  0: not
  significant, 1: significant and positive, -1: significant
  and negative}
}
\description{
  Identify any significant differences between predicted
  data before an impact event and predicted data after an
  impact event
}
\details{
  This function finds the differences for every predicted
  grid cell for every bootstrap replicate.  Quantiles are
  used to determine whether each difference is
  significantly different from zero and if so, in what
  direction.
}
\examples{
\dontrun{
getDifferences(beforePreds, afterPreds)}
}

