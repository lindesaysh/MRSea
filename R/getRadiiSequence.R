#' Function for obtaining a sequence of range parameters for the bivariate CReSS smoother
#' 
#' @param method One of "original" or "variogram". The original method is described in Scott-Hayward et al 2013 (reference below) and the variogram method is in the [vignette for variogram radii selection](https://lindesaysh.github.io/MRSea/articles/web/UserRadiiChoice_MRSea.html)
#' @param numberofradii The number of range parameters for SALSA to use when fitting the CReSS smooth.  The default is 8.  Remember, the more parameters the longer SALSA will take to find a suitable one for each knot location.
#' @param distMatrix  Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances. Euclidean distances created using \code{\link{makeDists}}.
#' @param basis character stating whether a 'gaussian' or 'exponential' basis is being used. 
#' #' @param rvin Optional if `method = "original"`. Two parameter vector stating the minimum and maximum range of r for a gaussian basis.
#' @param xydata Required if `method = "variogram". Data frame containing columns for x and y coordinates. x is assumed to be the first of the two columns
#' @param response Required if `method = "variogram"`. Vector of response values for use in `gstat::variogram`.  These values should be approximately normally distributed.
#' @param alpha Optional if `method = "variogram"`. Numeric parameter for the `gstat::variogram` function giving the direction in plane(x,y)
#' @param showplots (`default = FALSE`). Optional if `method = "variogram"`. If `TRUE` the output of `gstat::variogram` and `gstat::fit.variogram` are shown.
#' @param ... Other parameters for the `gstat::variogram` function. 
#' 
#' @details
#' The range parameter determines the range of the influence of each knot.  Small numbers indicate local influence and large ones, global influence.
#' 
#' @references
#' Scott-Hayward, L.; M. Mackenzie, C.Donovan, C.Walker and E.Ashe.  Complex Region Spatial Smoother (CReSS). Journal of computational and Graphical Statistics. 2013. 
#' DOI: 10.1080/10618600.2012.762920
#' 
#' [vignette for variogram radii selection](https://lindesaysh.github.io/MRSea/articles/web/UserRadiiChoice_MRSea.html)
#' 
#' @return
#' 
#' * Original method: This function returns a vector containing a sequence of range parameters.
#' 
#' * Variogram method: This function returns a vector containing a sequence of range parameters.  If an even number of radii is requested, this is reduced by one to give an odd length sequence where the middle number was the best range parameter from the variogram. The outputs of the variogram model can be found in the attributes of the returned object under `vg.fit`. 
#' 
#' **Note** If the estimated range from the variogram model exceeds the maximum distance in the distance matrix, the sequence reverts to the original method. A warning will be printed and the method used is an attribute in the returned sequence ("Method").
#' 
#' @examples
#' 
#' # load data
#' data(ns.data.re)
#' rad.dat <- dplyr::filter(ns.data.re, impact==0, Year==9, MonthOfYear == 3)
#' # load knot grid data
#' data(knotgrid.ns)
#' 
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(rad.dat$x.pos, rad.dat$y.pos), na.omit(knotgrid.ns))

#' # choose sequence of radii
#' r_seq<-getRadiiSequence(method = "variogram",
#'                         numberofradii = 8, 
#'                         xydata = rad.dat[,c("x.pos", "y.pos")], 
#'                         response = log(rad.dat$birds +1 ), 
#'                         basis = "gaussian", 
#'                         distMatrix = distMats$dataDist)
#' 
#' r_seq
#' attr(r_seq, "vg.fit")
#' attr(r_seq, "Method")
#' 
#'  r_seq<-getRadiiSequence(method = "original",
#'                          numberofradii = 8, 
#'                          distMatrix = distMats$dataDist, 
#'                          basis="gaussian")
#' 
#' r_seq
#' attr(r_seq, "Method")
#' 
#' @export
#' 


getRadiiSequence<-function(method=NULL, numberofradii=10, 
                           distMatrix, basis,
                           rvin=NULL,
                           xydata, response, 
                           alpha = 0, vgmmodel = "Sph", 
                           showplots = FALSE, 
                           ...){
  
  if(is.null(distMatrix)) stop("**** No distance matrix provided.****")
  
  ifelse(method %in% c("original", "variogram"),1 , stop("*** method not one of original or variogram"))
  
  distMatrix[which(is.infinite(distMatrix), arr.ind = T)]<-NA
  
  # ~~~~~~~ variogram method ~~~~~~~~~ # 
  
  if(method == "variogram"){
    data <- data.frame(xydata, response)
    names(data) <- c("x", "y", "response")
    
    sp::coordinates(data) = ~x + y
    suppressWarnings({
      vg <- gstat::variogram(response ~ 1, data, alpha= alpha, ...)
      fit.vg <- try(gstat::fit.variogram(vg, model=gstat::vgm(vgmmodel)), silent = TRUE)
      
      if(class(fit.vg)[1]=="try-error"){
        fit.vg <- gstat::fit.variogram(vg, model=gstat::vgm(vgmmodel), fit.method = 6)
      }
    })
    #print(fit.vg)
    if(showplots == TRUE){
      par(mfrow=c(1,2))
      print(plot(vg))
      print(plot(fit.vg, cutoff=max(vg$dist)))
    }
    best.s <- fit.vg$range[2]
    
    
    if(max(distMatrix, na.rm = TRUE)<best.s){
      method = "original"
      print(paste0("Range (", 
                   round(best.s, 2), 
                   ") is greater than maximum distance in distMatrix (",
                   round(max(distMatrix, na.rm = TRUE), 3), 
                   ") so r_seq created using method = 'original'"
      ))
    }else{
      vg.gap <- vg$dist[2] - vg$dist[1]
      nr <- floor(numberofradii/2)
      l <- seq(best.s - (vg.gap * nr), best.s, length=nr)
      
      if(min(l)<0){
        l <- seq(min(abs(l)), best.s, length=nr)
      }
      
      u <- seq(best.s, best.s + (vg.gap * nr), length=nr)
      s_seq <- unique(c(l, u))
      
      
      if(basis=='gaussian'){
        
        r_seq <- 1/((sqrt(2) * s_seq))
        attr(r_seq, "Method") <- "Variogram"
        attr(r_seq, "vg.fit") <- fit.vg
      }
      
      if(basis=='exponential'){
        r_seq <- sqrt(s_seq)
        attr(r_seq, "Method") <- "Variogram"
        attr(r_seq, "vg.fit") <- fit.vg
      }
    } 
  }
  
  # ~~~~~~~~~ Original method ~~~~~~~~~~~~~~~~~
  # from CRESS paper
  
  if(method == "original"){
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
    }
    
    if(basis=='exponential'){
      numberofradii = numberofradii+2
      # establish smallest observation-knot distance
      rmin<- sqrt(max(distMatrix, na.rm=TRUE)/21)
      rmax<- sqrt(max(distMatrix, na.rm=TRUE)/3e-7)
      r_seq <- exp(seq(log(rmin), log(rmax), length=numberofradii))[-c(1,numberofradii)]
    }
    attr(r_seq, "Method") <- "Original"
    
  }
  return(r_seq)
  
}


getRadiiChoices.vario<-function(numberofradii=10, xydata, response, 
                                basis, alpha = 0, vgmmodel = "Sph", 
                                showplots = FALSE, distMatrix = NULL, ...){
  
  if(is.null(distMatrix)) stop("**** No distance matrix provided.****")
  
  
  
  return(r_seq) 
}
