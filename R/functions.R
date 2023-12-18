#' Function for obtaining a sequence of range parameters for the CReSS smoother
#' 
#' @description
#' `r lifecycle::badge("superseded")`
#' `getRadiiChoices()` has been superseded in favour of `getRadiiSequence()`
#' 
#' @param numberofradii The number of range parameters for SALSA to use when fitting the CReSS smooth.  The default is 8.  Remember, the more parameters the longer SALSA will take to find a suitable one for each knot location.
#' @param distMatrix  Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances. Euclidean distances created using \code{\link{makeDists}}.
#' @param basis character stating whether a 'gaussian' or 'exponential' basis is being used. 
#' @param rvin Two parameter vector stating the minimum and maximum range of r for a gaussian basis. 
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
#' r_seq<-getRadiiChoices(8, distMats$dataDist, basis="gaussian")
#' 
#' @export
#' 


getRadiiChoices<-function(numberofradii=10, distMatrix, basis, rvin=NULL){
  
  
  if(is.null(distMatrix)) stop("**** No distance matrix provided.****")
  
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
    attr(r_seq, "Method") <- "Original"
    return(r_seq)  
  }
  
  if(basis=='exponential'){
    numberofradii = numberofradii+2
    # establish smallest observation-knot distance
    rmin<- sqrt(max(distMatrix, na.rm=TRUE)/21)
    rmax<- sqrt(max(distMatrix, na.rm=TRUE)/3e-7)
    r_seq <- exp(seq(log(rmin), log(rmax), length=numberofradii))[-c(1,numberofradii)]
    attr(r_seq, "Method") <- "Original"
    return(r_seq)
  }
}


#-----------------------------------------------------------------------------
#' Function for obtaining a sequence of range parameters for the CReSS smoother
#' 
#' @description
#' `r lifecycle::badge("superseded")`
#' `getRadiiChoices.vario()` has been superseded in favour of `getRadiiSequence()`
#' 
#' @param numberofradii The number of range parameters for SALSA to use when fitting the CReSS smooth.  The default is 8.  Remember, the more parameters the longer SALSA will take to find a suitable one for each knot location.
#' @param xydata Data frame containing columns for x and y coordinates. x is assumed to be the first of the two columns
#' @param response vector of response values for use in `gstat::variogram`.  These values should be approximately normally distributed.
#' @param basis character stating whether a 'gaussian' or 'exponential' basis is being used. 
#' @param alpha numeric parameter for the `gstat::variogram` function giving the direction in plane(x,y)
#' @param showplots (`default = FALSE`). If `TRUE` the output of `gstat::variogram` and `gstat::fit.variogram` are shown.
#' @param distMatrix Matrix of distances between data locations and knot locations (n x k). May be Euclidean or geodesic distances. Euclidean distances created using \code{\link{makeDists}}. This is used as a check to ensure the estimated range parameter does not exceed the maximum distance on the surface.  If it does, then the original \code{\link{getRadiiChoices}} function is used and the `distMatrix` parameter is a requirement for this. 
#' @param ... Other parameters for the `gstat::variogram` function. 
#' 
#' 
#' @details
#' The range parameter determines the range of the influence of each knot.  Small numbers indicate local influence and large ones, global influence.
#' 
#' @return
#' This function returns a vector containing a sequence of range parameters.  If an even number of radii is requested, this is reduced by one to give an odd length sequence where the middle number was the best range parameter from the variogram. 
#' The outputs of the variogram model can be found in the attributes of the returned object under `vg.fit`.
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
#' r_seq<-getRadiiChoices.vario(8, xydata = rad.dat[,c("x.pos", "y.pos")], 
#'                              response = log(rad.dat$birds +1 ), 
#'                              basis = "gaussian", 
#'                              distMatrix = distMats$dataDist)
#' 
#' r_seq
#' attr(r_seq, "vg.fit")
#' @export
#' 


getRadiiChoices.vario<-function(numberofradii=10, xydata, response, 
                                basis, alpha = 0, vgmmodel = "Sph", 
                                showplots = FALSE, distMatrix = NULL, ...){
  
  if(is.null(distMatrix)) stop("**** No distance matrix provided.****")
  
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
  
  distMatrix[which(is.infinite(distMatrix), arr.ind = T)]<-NA
  if(max(distMatrix)>best.s){
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
  }else{
    r_seq <- getRadiiChoices(numberofradii, distMatrix, basis=basis)
    attr(r_seq, "vg.fit") <- paste0("Range (", 
                                    round(best.s, 2), 
                                    ") is greater than maximum distance in distMatrix (",
                                    round(max(distMatrix), 3), 
                                    ") so r_seq created by getRadiiChoices()"
    )
  }
  
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
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function makes two Euclidean distance matrices.  One for the distances between all spatial observations and all spatial knot locations.  The other, if specified, is the distances between knot locations. It is experimental owing to the addition of the creation of infinity-block distance matrices. The original functionality is unchanged.
#' 
#' @param datacoords Coordinates of the data locations. If a design matrix interaction is required, add a third column containing the factor level. 
#' @param knotcoords Coordinates of the legal knot locations. If a design matrix interaction is required, add a third column containing the factor level for each knot location.
#' @param knotmat (\code{default=TRUE}). Should a matrix of knot-knot distances be created
#' @param polys (\code{default=NULL}). If geodesic distances are to be calculated, provide a list of polygons defining exclusion areas. 
#' @param type (\code{default='A'}).  One of 'A' or 'B'.  'A' is used when the \code{knotcoords} are a subset of \code{datacoords} AND the attributes of \code{knotcoords} give the index of points from \code{datacoords} (as happens if \code{getKnotgrid()} is used).  'B' is used for prediction (when \code{datacoords} is a prediction grid and so \code{knotcoords} is not a subset) or when the knotgrid was not generated using \code{getKnotgrid()}.
#' @param plot.transition (\code{default=TRUE}). Logical stating whether to plot the transition matrix.  Useful to see if the boundaries are being obeyed.
#' @param grid.dim This is a vector of length two which specifies the dimesions of the grid used to create the transition matrix (default \code{c(100, 100)}.  If the transition matrix shows that boundaries are being ignored, the grid dimensions will need to increase. However, increasing the grid, whilst improving accuracy, also increases computational time. 
#' 
#' @details
#' The data-knot matrix is used in the CReSS basis and the knot-knot matrix is used in SALSA to determine where a nearest knot to `move' should be.
#' 
#' If three columns are provided for `datacoords` and `knotcoords` the matrix returned has infinity for distances between knots and data associated with differing factor levels. 
#'  
#' @examples
#' # load data
#' data(ns.data.re)
#' # load knot grid data
#' data(knotgrid.ns)
#'  
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(ns.data.re$x.pos, ns.data.re$y.pos), na.omit(knotgrid.ns))
#' 
#' # ~~~~~~~~~~~~~~~~~~
#' # Example with block-infinity distance matrix
#' data(nysted.analysisdata)
#' myknots <- selectFctrKnots(nysted.analysisdata[,c('x.pos', 'y.pos', 'impact')], nk=150)
#' 
#' dists <- makeDists(datacoords = nysted.analysisdata[,c('x.pos', 'y.pos', 'impact')], 
#'                    knotcoords = myknots, 
#'                    knotmat = TRUE)
#' @export
#' 

makeDists<-function(datacoords, knotcoords, knotmat=TRUE, 
                        polys=NULL, type='A', plot.transition=FALSE, grid.dim = c(100, 100)){
  
  stopifnot(ncol(knotcoords) == ncol(datacoords), (ncol(datacoords) | ncol(knotcoords)) < 4)
  
  if(ncol(knotcoords)==3){
    knotcoords.factor = knotcoords[,3]
    if(!is.factor(knotcoords.factor)) knotcoords.factor <- as.factor(knotcoords.factor)
    knotcoords = knotcoords[,1:2]
    datacoords.factor = datacoords[,3]
    if(!is.factor(datacoords.factor)) datacoords.factor <- as.factor(datacoords.factor)
    datacoords = datacoords[, 1:2]
    factorknots = TRUE
  }else{
    factorknots = FALSE
  }
  
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
      }
    } # end model
  } # end geodesic
  
  if(factorknots){
    kcontr <- model.matrix(~-1 + knotcoords.factor)
    
    if(knotmat == TRUE){
      # make square (kxk)
      kcontr.m <- kcontr %*% t(kcontr)
      # make 0's into infs and make diagonals 0's
      kcontr.inf <- ifelse(kcontr.m == 0, Inf, kcontr.m)
      diag(kcontr.inf) <- 0
      # change distance matrix
      knotDist <- knotDist + kcontr.inf
    }
    
    dcontr <- model.matrix(~-1 + datacoords.factor)
    dkcontr <- dcontr %*% t(kcontr)
    dkcontr.inf <- ifelse(dkcontr == 0, Inf, dkcontr)
    #diag(dkcontr.inf) <- 0
    
    d2k <- d2k + dkcontr.inf
    
  }
  
  
  if(knotmat==T){
    return(list(dataDist=d2k, knotDist = knotDist))
  }else{
    return(list(dataDist=d2k))
  }
  
}


#-----------------------------------------------------------------------------
#' Get Pearsons residuals for point process model
#' 
#' This function calculates the residuals for a point process model.
#' 
#' @param model A point process model fitted in a glm framework
#' @param type Type of residual required. Choices are "Pearson" (default) or "raw"


getPPresiduals<-function(model, type='Pearson'){
  data<-model$data
  lambda<-fitted(model)
  z<-data$response==1
  
  if(type=='Pearson'){
    indicator<-(lambda > .Machine$double.eps)
    rp.dens<-(-indicator) * sqrt(lambda)
    rp.disc<-1/(sqrt(lambda[z]))
  }
  if(type=='raw'){
    rp.dens <- (-lambda)
    rp.disc <- rep.int(1, sum(z))
  }
  nquad<-length(data$response==0)
  discretepad <- numeric(nquad)
  discretepad[z] <- rp.disc
  #wt <- w.quad(Q)
  val <- discretepad + ((model$model$'(weights)') * rp.dens)
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


#-----------------------------------------------------------------------------
#' Get block leave-one-out CV score for point process model
#' 
#' This function calculates a leave-one-out cross-validation score for point process models
#' 
#' @param model.obj a point process gamMRSea model
#' @param proj4string string giving CRS of x.pos and y.pos in date in `model.obj`
#' @param blocksize size of the blocks, give in m if `proj4string` is in degrees or m, or in km if data in  `model.obj`/`proj4string` is in km
#' @param folds (default = 10). number stating many folds for CV
#' @param replicates (default = 100). number stating how many times to replicate the CV score
#' @param showProgress (default = FALSE) logical stating whether to print progress bar

# getPPCV <- function(model.obj, proj4string, blocksize, 
#                     folds = 10, replicates = 100, showProgress = FALSE){
#   
#   if(class(model.obj)[1]!='gamMRSea')stop('Model object not of class gamMRSea')
#   if(model.obj$splineParams[[1]]$modelType!='pointProcess')stop('Model object not a point process model')
#   
#   
#   dat<-model.obj$data
#   spdata<-sp::SpatialPointsDataFrame(coords = cbind(dat$x.pos, dat$y.pos), 
#                                      data = dat,
#                                      proj4string = proj4string)
#   
#   #offset <- seq(0.1, 0.9, by = 0.1)
#   scores<-vector(length=replicates)
#   
#   for(r in 1:replicates){
#     if(showProgress) {if((r/10)%%1 == 0){cat(b, '\n')}else{cat('.')}}
#     
#     test<-try(blockCV::spatialBlock(spdata, species='response', k=folds, 
#                                     theRange = blocksize,
#                                 xOffset = 0.5, #sample(offset,1),
#                                 yOffset = 0.5, #sample(offset, 1),
#                                 verbose=FALSE, progress=FALSE, 
#                                 showBlocks = FALSE),
#               silent = TRUE)
#     
#     if(class(test)=='try-error'){
#       scores[r]<-NA
#     }else{
#       # match folds to data
#       dat<-dat %>%
#         mutate(cvfolds = test$foldID,
#                index = 1:nrow(dat))
#       
#       # vector of folds
#       cv.folds <- sort(unique(dat$cvfolds))
#       # make object to store squared resid
#       lo_rss <- vector(length=length(cv.folds))
#       
#       for(b in cv.folds){
#         
#         # # select training data
#         # NOT NEEDED, CHANGE WEIGHTS NOT DATASET SIZE
#         # lo_dat0 <- filter(dat, response==0)
#         # lo_dat1 <- filter(dat, response>=1, cvfolds!=b)
#         # lo_dat <- rbind(lo_dat1, lo_dat0)
#         
#         cvdat<-dat
#         cvdat$pp.wts.f <- cvdat$pp.wts
#         cvdat$pp.wts.f[dat$cvfolds==b] <- 0
#         
#         # ggplot() + 
#         #   geom_tile(data=filter(dat, response==0), 
#         #             aes(x.pos, y.pos, fill=as.factor(cvblocks)),
#         #             height=2, width=2, alpha=1/2) + 
#         #   coord_equal() + theme_bw() + 
#         #   geom_point(data=filter(dat, response==1), 
#         #             aes(x.pos, y.pos, colour=as.factor(cvblocks)),
#         #            size=1, alpha=1)
#         # 
#         # ggplot() + 
#         #   geom_tile(data=filter(dat, response==0), 
#         #             aes(x.pos, y.pos, fill=as.factor(cvblocks)),
#         #             height=2, width=2, alpha=1/2, show.legend = FALSE) + 
#         #   coord_equal() + theme_bw() + 
#         #   geom_point(data=filter(dat, response==1), 
#         #              aes(x.pos, y.pos),
#         #              size=1, alpha=1) + 
#         #   geom_point(data=filter(lo_dat, response==1), 
#         #              aes(x.pos, y.pos, colour=as.factor(cvblocks)),
#         #              size=1, alpha=1, show.legend = FALSE)
#         
#         # # update distance matrix for reduced dataset
#         #NOT NEEDED AS DATASIZE NOT CHANGING
#         # sp<-model.obj$splineParams
#         # sp[[1]]$dist <- model.obj$splineParams[[1]]$dist[lo_dat$index,]
#         # 
#         # update model for reduced data set
#         lom_model <- update.gamMRSea(model.obj, .~., data=cvdat,weights=pp.wts.f)
#         
#         # get predicted intensity for all points
#         mydat_temp = dat %>% 
#           mutate(preds =  predict.gamMRSea(object = lom_model, newdata = dat, g2k = model.obj$splineParams[[1]]$dist))
#         
#         # for each cv fold, find sum of presences (observed count)
#         cv.fold.intensity<-mydat_temp %>%
#           group_by(cvfolds) %>% 
#           summarise(npts = sum(response))
#         
#         # for each cvfold, find the sum of the estimated intensity for quad points
#         cv.fold.quads.intensity<-filter(mydat_temp, response==0) %>%
#           group_by(cvfolds) %>% 
#           summarise(npts = sum(preds))
#         
#         # join the two together to get data frame of model fit and observed count for
#         # each fold and calculate absolute residual. 
#         # select only the validation fold b
#         modelfit<-left_join(cv.fold.intensity, cv.fold.quads.intensity, by="cvfolds") %>%
#           mutate(resids = abs(npts.x - npts.y)) %>%
#           filter(cvfolds == b)
#         
#         # find squared residual and append for all b
#         lo_rss[b] <- modelfit$resids**2
#       } # end k loop
#       
#       scores[r]<-mean(lo_rss, na.rm=TRUE)
#     }
#     
#   }
#   return(list(cv = mean(scores, na.rm=TRUE), scores = scores))
# }



getCVpp<-function(model.obj, folddata, kfold=10, replicates=1){
  
  if(class(model.obj)[1]!='gamMRSea')stop('Model object not of class gamMRSea')
  if(model.obj$splineParams[[1]]$modelType!='pointProcess')stop('Model object not a point process model')
  
  scores<-matrix(NA, nrow=replicates, ncol=kfold)
  
  for(r in 1:replicates){
    if(replicates>1){
      myfold <- folddata[[r]]
    }else{
      myfold <- folddata
    }
    
    # make object to store squared resid
    sqresid = absresid = vector(length=length(kfold))
    
    for(f in 1:kfold){
      
      dat <- model.obj$data
      # select training data
      # set weights for fold data to zero to be equivalent to thinning
      dat$pp.wts.f <- dat$pp.wts
      dat$pp.wts.f[myfold[[f]][[2]]] <- 0
      # 
      # update model for reduced data set
      lom_model <- update.gamMRSea(model.obj, .~., data=dat,weights=pp.wts.f)
      
      # get predicted intensity for all points
      mydat_temp = dat %>% 
        mutate(preds =  predict.gamMRSea(object = lom_model, newdata = dat, g2k = model.obj$splineParams[[1]]$dist))
      
      # for each cv fold, find sum of presences (observed count)
      cv.fold.intensity<-sum(mydat_temp$response[myfold[[f]][[2]]])
      
      # for each cvfold, find the sum of the estimated intensity for quad points
      cv.fold.quads.intensity<-sum(mydat_temp$preds[which(mydat_temp$response == 0 &  mydat_temp$pp.wts.f==0)])
      
      # join the two together to get data frame of model fit and observed count for
      # each fold and calculate absolute residual. 
      # select only the validation fold b
      resid = cv.fold.intensity - cv.fold.quads.intensity
      
      # find squared residual and append for all b
      sqresid[f] <- resid**2
      #absresid[f] <- abs(resid)
    } # end k loop
    
    scores[r,]<-sqresid
    
  }
  
  allcvs <- apply(scores, 1, mean)
  cv = mean(allcvs)
  
  if(replicates>1){
    cis<-quantile(allcvs, probs = c(0.025,0.975))
    return(list(cv = cv, cis = cis, scores = scores))
  }else{
    return(list(cv = cv, scores = scores))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#-----------------------------------------------------------------------------
#' Select candidate knots for multi-level factor interaction
#' 
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' This function finds a number of starting knot locations for different factor levels and is used for creating an interaction term. 
#' 
#' @param data Coordinates of the data locations in columns 1 and 2 and the associated factor level in the third. 
#' @param nk  The number of knots to be selected for each factor level.  This may be a single number (for all factor levels) or a vector of numbers to allow a different number of knots to be selected for each factor level.
#' @param s.eed `default = 1`. Set the seed for the selection to ensure reproducibility. 
#' 
#' @details
#' The function returns a data frame with three columns; x and y locations for candidate knots and associated factor level. 
#' 
#' @examples
#' # load data
#' data(nysted.analysisdata)
#' 
#' myknots <- selectFctrKnots(nysted.analysisdata[,c("x.pos", "y.pos", "impact")], 
#'                            nk=150)
#' 
#' @export
#' 

selectFctrKnots <- function(data, nk, s.eed = 1){
  
  fctrs = unique(data[,3])
  n.level = length(fctrs)
  myknots = NULL
  
  if(length(nk)==1){
    nk <- rep(nk, n.level)
  }
  
  for(i in 1:n.level){
    faclevel<-fctrs[i]
    set.seed(s.eed)
    tempkg<-getKnotgrid(coordData = data[data[,3] == faclevel, 1:2], 
                        numKnots = nk[i], plot = FALSE)
    myknots<-rbind(myknots, data.frame(tempkg, faclevel = faclevel))
  }
  
  names(myknots) <- names(data)
  row.names(myknots) <- NULL
  
  return(myknots)
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#-----------------------------------------------------------------------------
#' Select starting knots for multi-level factor interaction
#' 
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' This function finds a number of starting knot locations for different factor levels and is used for creating an interaction term. 
#' 
#' @param knots Coordinates of the knot locations in columns 1 and 2 and the associated factor level in the third. 
#' @param nk  The number of knots to be selected for each factor level.  This may be a single number (for all factor levels) or a vector of numbers to allow a different number of knots to be selected for each factor level.
#' @param s.eed `default = 1`. Set the seed for the selection to ensure reproducibility. 
#' 
#' @details
#' The function returns the row number of the `knots` parameter for the selected knots. 
#'
#' @examples
#' 
#'#' # load data
#' data(nysted.analysisdata)
#' 
#' myknots <- selectFctrKnots(nysted.analysisdata[,c('x.pos', 'y.pos', 'impact')], 
#'                            nk=150)
#' startknotlocs <- selectFctrStartk(myknots, nk=5, s.eed = 4)
#' 
#' 
#' @export
#' 

selectFctrStartk <- function(knots, nk, s.eed = 1){
  
  n.level = length(unique(knots[,3]))
  startknotlocs = NULL
  
  if(length(nk)==1){
    nk <- rep(nk, n.level)
  }
  
  for(i in 1:n.level){
    dat <- knots[knots[,3] == unique(knots[,3])[i], 1:2]
    set.seed(s.eed)
    ks <- as.numeric(rownames(getKnotgrid(dat, 
                                          numKnots = nk[i], 
                                          plot = FALSE)))
    #ks <- ks + ((i-1)* nrow(dat))
    startknotlocs <- c(startknotlocs, ks)
  }
  return(startknotlocs)
} 



