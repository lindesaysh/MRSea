\name{return.reg.spline.fit.2d}
\alias{return.reg.spline.fit.2d}
\title{Code for adaptively spacing knots for a spatial smooth. The smoothing process uses a CReSS basis.}
\usage{
  return.reg.spline.fit.2d(splineParams, startKnots,
    winHalfWidth, fitnessMeasure = "BIC",
    maxIterations = 100, tol = 0, baseModel = NULL,
    radiusIndices = NULL, initialise = TRUE,
    initialKnots = NULL, interactionTerm = NULL,
    knot.seed = 10)
}
\arguments{
  \item{splineParams}{List object where the first element
  [[1]] contains a list of objects for the 2D SALSA fitting
  process: \code{knotDist}, \code{radii}, \code{dist},
  \code{gridresp}, \code{grid}, \code{datacoords},
  \code{response}, \code{knotgrid}, \code{minKnots},
  \code{maxKnots}, \code{gap}}

  \item{startKnots}{number of space-filled knots to start
  with (between minKnots and maxKnots)}

  \item{fitnessMeasure}{(default=BIC). Measure used to
  evaluate the fit. Other options are AIC, AICc, BIC, QICb
  (Quasi-Likelihood Information Criterion with log(n)
  penalty)}

  \item{maxIterations}{exchange/improve heuristic will
  terminate after \code{maxIterations} if still running}

  \item{baseModel}{starting model for SALSA to use.  Must
  not already contain a spatial smooth.}

  \item{radiusIndices}{vector of length \code{startKnots}
  identifying which radii (\code{splineParams[[1]]$radii})
  will be used to initialise the model}

  \item{initialise}{(default = TRUE). Logical stating
  whether or not to start with space-filled knots (TRUE) or
  user specified locations (FALSE)}

  \item{initialKnots}{If \code{initialise=FALSE} then the
  start locations for the knots are specified in
  \code{initialKnots}. Must be coordinates.}

  \item{interactionTerm}{(default=NULL). Specifies which
  term in \code{baseModel} the spatial smooth will interact
  with.  If NULL no interaction term is fitted}

  \item{winHalfWidth}{Half-width of window used to
  calculate region with biggest average residual magnitude}

  \item{tol}{Tolerance for difference between fit measures.
  E.g. \code{tol=2} means that the calculated fitness
  measures must be 2 units apart to be considered
  different}
}
\description{
  Code for adaptively spacing knots for a spatial smooth.
  The smoothing process uses a CReSS basis.
}
\details{
  The following are the details of the splineParams[[1]]
  objects.  Note.  If salsa1D has been run then details for
  those covariates will sit in splineParams[[2]] and
  onward.

  \code{knotDist}: matrix of knot to knot distances (k x
  k).  May be Euclidean or geodesic distances. Must be
  square and the same dimensions as
  \code{nrows(na.omit(knotgrid))}

  \code{radii} Sequence of range parameters for the CReSS
  basis from local (small) to global (large).  Determines
  the range of the influence of each knot.

  \code{dist}: matrix of distances between data locations
  and knot locations (n x k). May be Euclidean or geodesic
  distances.

  \code{gridresp} The first column of knotgrid

  \code{grid} Index of knotgrid locations.  Should be same
  length as \code{knotgrid} but with x=integer values from
  1 to number of unique x-locations and y= integer values
  from 1 to number of unique y-locations.

  \code{datacoords}: Coordinates of the data locations

  \code{response}: vector of response data for the
  modelling process

  \code{knotgrid}: grid of legal knot locations.  Must be a
  regular grid with c(NA, NA) for rows with an illegal knot

  \code{minKnots}: minimum number of knots to fit

  \code{maxKnots}: maximum number of knots to fit

  \code{gap}: Minimum gap between knots (in unit of
  measurement of \code{datacoords})
}

