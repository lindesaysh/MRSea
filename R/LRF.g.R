#' Function for creating an Gaussian basis function for a spatial smooth using the CReSS method.
#' 
#' This function calculates a local radial Gausiian basis matrix for use in \code{\link{runSALSA2D}}.
#' 
#' @param radiusIndices Vector of length startKnots identifying which radii (splineParams[[1]]$radii) will be used to initialise the model
#' @param dists Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances.
#' @param radii Sequence of range parameters for the CReSS basis from local (small) to global (large). Determines the range of the influence of each knot.
#' @param aR Index of knot locations. The index contains numbers selected by SALSA from 1 to the number of legal knot locations \code{na.omit(knotgrid)}. Used to specify which columns of \code{dists} should be used to construct the basis matrix.
#' 
#' @details
#' Calculate a local radial basis matrix for use in \code{\link{runSALSA2D}}.  The distance matrix input may be Euclidean or geodesic distances.
#' 
#' @return
#' Returns a basis matrix with one column for each knot in \code{aR} and one row for every observation (i.e. same number of rows as \code{dists})
#' 
#' @examples
#' 
#' # load data
#' data(ns.data.re)
#' # load knot grid data
#' data(knotgrid.ns)
#' 
#' splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour'))
#' 
#' #set some input info for SALSA
#' ns.data.re$response<- ns.data.re$birds
#' 
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns), knotmat=FALSE)
#' 
#' # choose sequence of radii
#' r_seq<-getRadiiChoices(8, distMats$dataDist)
#'
#' # using the fourth radius and picking 5 knots
#' basis<-LRF.g(radiusIndices=rep(4, 5), dists=distMats$dataDist, radii = r_seq, 
#'         aR=c(3, 10, 15, 28, 31))
#' 
#' @export
#' 
LRF.g<- function(radiusIndices, dists, radii,aR){
  
  
  for (i in 1:length(aR)){
    zhold<- dists[,aR[i]]
    r<- radii[radiusIndices[i]]
  
    zhold<- exp(-(zhold*radii[radiusIndices[i]])**2)   
  
    if (i==1) {B<- zhold} else {B<- cbind(B,zhold)}
  }
  B <- data.frame(B)
  for(i in 1:ncol(B)){
    names(B)[i] <- paste('b', i, sep='')
  }
  B=as.matrix(B)
  #B<-scale(B, center=T, scale=F)
  return(B)
}