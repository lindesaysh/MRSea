#-----------------------------------------------------------------------------
#' find the plotting dimensions for quilt.plot when using a regular grid
#' 
#' @param x.pos Vector of x-coordinates in dataset
#' @param y.pos Vector of y-coordinates in dataset
#' @param segmentWidth Width of each grid cell of data (in same units as \code{x.pos})
#' @param segmentLength Length of each grid cell of data (in same units as \code{y.pos})
#' 
#' 
#' @examples getPlotdimensions(data$x.pos, data$y.pos, segmentWidth=500, segmentLength=500)
#' 
#' @export
#' 
getPlotdimensions<-function(x.pos, y.pos, segmentWidth, segmentLength){
  nrow<-(range(x.pos)[2] - range(x.pos)[1])/segmentWidth
  ncol<-(range(y.pos)[2] - range(y.pos)[1])/segmentLength
  return(c(nrow, ncol))
}