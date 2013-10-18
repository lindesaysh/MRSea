\name{getPlotdimensions}
\alias{getPlotdimensions}
\title{find the plotting dimensions for quilt.plot when using a regular grid}
\usage{
  getPlotdimensions(x.pos, y.pos, segmentWidth,
    segmentLength)
}
\arguments{
  \item{x.pos}{Vector of x-coordinates in dataset}

  \item{y.pos}{Vector of y-coordinates in dataset}

  \item{segmentWidth}{Width of each grid cell of data (in
  same units as \code{x.pos})}

  \item{segmentLength}{Length of each grid cell of data (in
  same units as \code{y.pos})}
}
\description{
  find the plotting dimensions for quilt.plot when using a
  regular grid
}
\examples{
getPlotdimensions(data$x.pos, data$y.pos, segmentWidth=500, segmentLength=500)
}

