
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
# getRadiiChoices<-function(numberofradii=8, distMatrix){
#   numberofradii = numberofradii+2
#   # establish smallest observation-knot distance
#   rmin<- sqrt(max(distMatrix)/21)
#   rmax<- sqrt(max(distMatrix)/3e-7)
#   r_seq <- exp(seq(log(rmin), log(rmax), length=numberofradii))[-c(1,numberofradii)]
#   return(r_seq)
# }

# getRadiiChoices<-function(numberofradii=10, distMatrix){
#   #mindist<-mean(apply(distMatrix, 2, min))
#   #maxdist<-mean(apply(distMatrix, 2, max))
#   #r_seq=1/(seq((mindist), (maxdist),length=numberofradii+1))
#   mindist<-mean(apply(distMatrix, 2, min))
#   maxdist<-mean(apply(distMatrix, 2, max))
#   r_seq=exp(seq(log(1/mindist)*1.5, log(1/maxdist),length=10))
#   return(r_seq[-1])
# }

# getRadiiChoices<-function(numberofradii=10, distMatrix){
#   
#   mindist<-mean(apply(distMatrix, 2, min))
#   maxdist<-mean(apply(distMatrix, 2, max))
#   r_seq=exp(seq(log(1/mindist), log(1/maxdist)*1.5,length=20))
#   
#   aRsel<- sample(1:ncol(distMatrix), 5)
#   r_min=r_max=vector(length = length(aRsel))
#   for(a in 1:5){
#     b<-LocalRadialFunction(radiusIndices = c(1:length(r_seq)),dists = distMatrix, radii = r_seq, aR = rep(aRsel[a], length=length(r_seq)))
#     means<-apply(b, 2, mean)
#     r_min[a]<-r_seq[(which(means>0.01)[1])]
#     r_max[a]<-r_seq[(which(means>0.90)[1])]
#     
#   }
#   r_min<-mean(r_min)
#   r_max<-mean(r_max)
#   r_seq=exp(seq(log(r_min), log(r_max),length=numberofradii))
#   #b<-LocalRadialFunction(radiusIndices = c(1:length(r_seq)),dists = distMatrix, radii = r_seq, aR = rep(aRsel[a], length=length(r_seq)))
# #   summary(b)
# #   par(mfrow=c(2,2))
# #   for(i in 1:(length(r_seq))){
# #     print(i)
# #     quilt.plot(data$x.pos, data$y.pos, b[,i], asp=1, zlim=c(0,1))
# #   }
#   
#   return(r_seq)
# }

getRadiiChoices<-function(numberofradii=10, distMatrix){
  
  minDist <- mean(apply(distMatrix,2,min))
  meanDist <- mean(apply(distMatrix,2,mean))
  rval_max<- sqrt(-log(0.7)/meanDist**2)
  rval_min<- sqrt(-log(0.001)/meanDist**2)
  r_seq<- exp(seq(log(rval_min), log(rval_max), length=numberofradii))
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

