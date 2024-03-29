% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/return.reg.spline.fit.2d.R
\name{return.reg.spline.fit.2d}
\alias{return.reg.spline.fit.2d}
\title{Wrapper function for running SALSA2D}
\usage{
return.reg.spline.fit.2d(
  splineParams,
  startKnots,
  winHalfWidth,
  fitnessMeasure = "BIC",
  maxIterations = 10,
  tol = 0,
  baseModel = NULL,
  radiusIndices = NULL,
  initialise = TRUE,
  initialKnots = NULL,
  initialaR = NULL,
  interactionTerm = NULL,
  knot.seed = 10,
  plot = FALSE,
  cv.opts,
  basis
)
}
\description{
Wrapper function for running SALSA2D
}
\author{
Cameron Walker, Department of Enginering Science, University of Auckland.
}
