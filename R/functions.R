#' Function for obtaining a sequence of range parameters for the CReSS smoother
#' 
#' @param numberofradii The number of range parameters for SALSA to use when fitting the CReSS smooth.  The default is 8.  Remember, the more parameters the longer SALSA will take to find a suitable one for each knot location.
#' @param distMatrix  Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances. Euclidean distances created using \code{\link{makeDists}}.
#' @param basis character stating whether a 'gaussian' or 'exponential' basis is being used. 
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

getRadiiChoices<-function(numberofradii=10, distMatrix, basis, rvin=NULL){
  
  distMatrix[which(is.infinite(distMatrix), arr.ind = T)]<-NA
  
  if(basis=='gaussian'){
    minDist <- mean(apply(distMatrix,2,min, na.rm=TRUE))
    meanDist <- mean(apply(distMatrix,2,mean, na.rm=TRUE))
    if (is.null(rvin)) {
      rval_max<- sqrt(-log(0.7)/meanDist**2)
      rval_min<- sqrt(-log(0.001)/meanDist**2)
    } else {
      rval_max = rvin[1]
      rval_min <- rvin[2]
    }
    r_seq<- exp(seq(log(rval_min), log(rval_max), length=numberofradii))
    return(r_seq)  
  }
  
  if(basis=='exponential'){
    numberofradii = numberofradii+2
      # establish smallest observation-knot distance
      rmin<- sqrt(max(distMatrix, na.rm=TRUE)/21)
      rmax<- sqrt(max(distMatrix, na.rm=TRUE)/3e-7)
      r_seq <- exp(seq(log(rmin), log(rmax), length=numberofradii))[-c(1,numberofradii)]
      return(r_seq)
  }
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
#' @param polys (\code{default=NULL}). If geodesic distances are to be calculated, provide a list of polygons defining exclusion areas. 
#' @param type (\code{default='A'}).  One of 'A' or 'B'.  'A' is used when the \code{knotcoords} are a subset of \code{datacoords} AND the attributes of \code{knotcoords} give the index of points from \code{datacoords} (as happens if \code{getKnotgrid()} is used).  'B' is used for prediction (when \code{datacoords} is a prediction grid and so \code{knotcoords} is not a subset) or when the knotgrid was not generated using \code{getKnotgrid()}.
#' @param plot.transition (\code{default=TRUE}). Logical stating whether to plot the transition matrix.  Useful to see if the boundaries are being obeyed.
#' @param grid.dim This is a vector of length two which specifies the dimesions of the grid used to create the transition matrix (default \code{c(100, 100)}.  If the transition matrix shows that boundaries are being ignored, the grid dimensions will need to increase. However, increasing the grid, whilst improving accuracy, also increases computational time. 
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

makeDists<-function(datacoords, knotcoords, knotmat=TRUE, polys=NULL, type='A', plot.transition=FALSE, grid.dim = c(100, 100)){
  
 if(is.null(polys)){
   # Euclidean
  if(length(which(is.na(knotcoords[,1])))>0) stop('remove NAs from knotcoords')
  
  d2k<- matrix(0, ncol=dim((knotcoords))[1], nrow=length(datacoords[,1]))
  for(i in 1:dim(knotcoords)[1]){
    d2k[,i]<- sqrt((datacoords[,1]-knotcoords[i,1])**2 + (datacoords[,2]-knotcoords[i,2])**2)
  }
  
  if(knotmat==T){
    #specify the knot-to-knot distances (this cannot have any NAs included); size k x k
    knotDist = as.matrix(dist(na.omit(knotcoords), method = "euclidean", diag = TRUE, upper=TRUE))
    return(list(dataDist=d2k, knotDist = knotDist))
  }else{
    return(list(dataDist=d2k))
  }
} # end euclidean
  
  # Geodesic
  if(!is.null(polys)){
    
    Nx=seq(min(c(datacoords[,1], knotcoords[,1])), max(c(datacoords[,1], knotcoords[,1])), length=grid.dim[1])
    Ny=seq(min(c(datacoords[,2], knotcoords[,2])), max(c(datacoords[,2], knotcoords[,2])), length=grid.dim[2])
  
    xygrid<-expand.grid(x=Nx, y=Ny)
  
    if(type=='B'){
      if(!is.null(names(datacoords))){
        knotcoords<-data.frame(knotcoords)
        names(knotcoords)<-names(datacoords)
      }
      datacoords2 = rbind(datacoords, knotcoords)
      geodistsoutput<-getGeoDist(xygrid=xygrid, polys=polys, datalocations=datacoords2, plot.transition=plot.transition)
      # select out knot columnns to get d2k
      d2k<-geodistsoutput$distance[(1:nrow(datacoords)),((nrow(datacoords)+1):nrow(datacoords2))]
      if(knotmat==T){
        #specify the knot-to-knot distances (this cannot have any NAs included); size k x k
        knotDist = geodistsoutput$distance[(nrow(datacoords)+1):nrow(datacoords2),(nrow(datacoords)+1):nrow(datacoords2)] 
        return(list(dataDist=d2k, knotDist = knotDist))
      }else{
        return(list(dataDist=d2k))
      }
    } # end prediction
    if(type=='A'){
      # data to data matrix
      geodistsoutput<-getGeoDist(xygrid=xygrid, polys=polys, datalocations=datacoords, plot.transition=plot.transition)
      # select out knot columnns to get d2k
      d2k<-geodistsoutput$distance[,attr(knotcoords, 'points.selected')]
    if(knotmat==T){
      #specify the knot-to-knot distances (this cannot have any NAs included); size k x k
      knotDist = geodistsoutput$distance[attr(knotcoords, 'points.selected'),attr(knotcoords, 'points.selected')] 
      return(list(dataDist=d2k, knotDist = knotDist))
    }else{
      return(list(dataDist=d2k))
    }
    } # end model
  } # end geodesic
}


#-----------------------------------------------------------------------------
#' Get Pearsons residuals for point process model
#' 
#' This function calculates the residuals for a point process model.
#' 
#' @param model A point process model fitted in a glm framework


getPPresiduals<-function(model){
  data<-model$data
  lambda<-fitted(model)
  z<-data$response==1
  indicator<-(lambda > .Machine$double.eps)
  rp.dens<-(-indicator) * sqrt(lambda)
  rp.disc<-1/(sqrt(lambda[z]))
  
  nquad<-length(data$response==0)
  discretepad <- numeric(nquad)
  discretepad[z] <- rp.disc
  #wt <- w.quad(Q)
  val <- discretepad + data$pp.wts * rp.dens
  return(val)
}


#-----------------------------------------------------------------------------
#' Get index of n largest residuals
#' 
#' This function evaluates the 5 largest residuals from a model.
#' 
#' @param model a gamMRSea model
#' @param n number of largest residuals to return

getlargestresid<-function(model, n=5){
  if (isS4(model)){
    indexdat<-order(rowSums(abs(residuals(model, type='pearson'))), decreasing = TRUE)[1:n]
  }else{
    data<-model$data
    if(model$splineParams[[1]]$modelType=='pointProcess'){
      val<-getPPresiduals(model)
      indexdat<-order(abs(val), decreasing = TRUE)[1:n]
    }else {
      indexdat<-order(abs(residuals(model, type='pearson')), decreasing = TRUE)[1:n]
    }  
  }
  return(indexdat)
}