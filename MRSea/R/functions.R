
#-----------------------------------------------------------------------------
#' Function for obtaining a sequence of range parameters for the CReSS smoother
#' 
#' @param numberofradii The number of range parameters for SALSA to use when fitting the CReSS smooth.  The default is 8.  Remember, the more parameters the longer SALSA will take to find a suitable one for each knot location.
#' @param distMatrix  Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances. Euclidean distances created using \code{\link{makeDists}}.
#' 
#' @details
#' The range parameter determines the range of the influence of each knot.  Small numbers indicate local influence and large ones, global influence.  
#' 
#' @references
#' Scott-Hayward, L.; M. Mackenzie, C.Donovan, C.Walker and E.Ashe.  Complex Region Spatial Smoother (CReSS). Journal of computational and Graphical Statistics. 2013. 
#' DOI: 10.1080/10618600.2012.762920
#' 
#' @return
#' This function returns a vector containing a sequence of range parameters.
#' 
#' @examples
#' 
#' # load data
#' data(ns.data.re)
#' # load knot grid data
#' data(knotgrid.ns)
#' 
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))
#'
#' # choose sequence of radii
#' r_seq<-getRadiiChoices(8, distMats$dataDist)
#' 
#' @export
#' 
getRadiiChoices<-function(numberofradii=8, distMatrix){
  numberofradii = numberofradii+2
  # establish smallest observation-knot distance
  rmin<- sqrt(max(distMatrix)/21)
  rmax<- sqrt(max(distMatrix)/3e-7)
  r_seq <- exp(seq(log(rmin), log(rmax), length=numberofradii))[-c(1,numberofradii)]
  return(r_seq)
}

#-----------------------------------------------------------------------------
#' Factor level response check
#'  
#' This function checks that there are some non-zero counts in each level of each factor variable for consideration in a model
#' 
#' @param factorlist Vector of factor variables specified in \code{model}.  Specified so that a check can be made that there are non-zero counts in all levels of each factor. 
#' @param data Data frame containing columns of covariates listed in \code{factorlist}.  Column names must match with names in \code{factorlist}
#' @param response A vector of response values
#' 
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' 
#' checkfactorlevelcounts(factorlist=c('floodebb', 'impact'), ns.data.re, 
#'      ns.data.re$birds)
#' 
#' @export
#' 
checkfactorlevelcounts<-function(factorlist, data, response){
  for(i in 1:length(factorlist)){
    
    pos<- which(names(data)==factorlist[i])
    minCount<- min(tapply(response,data[,pos], sum ))
    
    if(minCount<1){stop(paste(names(data)[pos], " cannot be fitted as a factor variable, due to zero counts for some levels", sep=""))}
    
    if(minCount>1){print(paste(names(data)[pos], " will be fitted as a factor variable; there are non-zero counts for all levels", sep=""))}  
  }
}

#-----------------------------------------------------------------------------
#' Make Euclidean distance matrices for use in CReSS and SALSA model frameworks
#' 
#' This function makes two Euclidean distance matrices.  One for the distances between all spatial observations and all spatial knot locations.  The other, if specified, is the distances between knot locations.
#' 
#' @param datacoords Coordinates of the data locations
#' @param knotcoords Coordinates of the legal knot locations
#' @param knotmat (\code{default=TRUE}). Should a matrix of knot-knot distances be created
#' 
#' @details
#' The data-knot matrix is used in the CReSS basis and the knot-knot matrix is used in SALSA to determine where a nearest knot to `move' should be.
#'  
#' @examples
#' # load data
#' data(ns.data.re)
#' # load knot grid data
#' data(knotgrid.ns)
#'  
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))
#' @export
#' 

makeDists<-function(datacoords, knotcoords, knotmat=TRUE){
  
  if(length(which(is.na(knotcoords[,1])))>0) stop('remove NAs from knotcoords')
  
  DistMatrixNaive<- matrix(0, ncol=dim((knotcoords))[1], nrow=length(datacoords[,1]))
  for(i in 1:dim(knotcoords)[1]){
    DistMatrixNaive[,i]<- sqrt((datacoords[,1]-knotcoords[i,1])**2 + (datacoords[,2]-knotcoords[i,2])**2)
  }
  if(knotmat==T){
    #specify the knot-to-knot distances (this cannot have any NAs included); size k x k
    knotDist = as.matrix(dist(na.omit(knotcoords), method = "euclidean", diag = TRUE, upper=TRUE))
    return(list(dataDist=DistMatrixNaive, knotDist = knotDist))
  }else{
    return(list(dataDist=DistMatrixNaive))
  }
}

