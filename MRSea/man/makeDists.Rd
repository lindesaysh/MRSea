\name{makeDists}
\alias{makeDists}
\title{Make Euclidean distance matrices for use in CReSS and SALSA model frameworks}
\usage{
  makeDists(datacoords, knotcoords, knotmat = TRUE)
}
\arguments{
  \item{datacoords}{Coordinates of the data locations}

  \item{knotcoords}{Coordinates of the legal knot
  locations}

  \item{knotmat}{(\code{default=TRUE}). Should a matrix of
  knot-knot distances be created}
}
\description{
  This function makes two Euclidean distance matrices.  One
  for the distances between all spatial observations and
  all spatial knot locations.  The other, if specified, is
  the distances between knot locations.
}
\details{
  The data-knot matrix is used in the CReSS basis and the
  knot-knot matrix is used in SALSA to determine where a
  nearest knot to `move' should be.
}
\examples{
# load data
data(ns.data.re)
# load knot grid data
data(knotgrid.ns)

# make distance matrices for datatoknots and knottoknots
distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))
}

