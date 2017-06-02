#' Bootstrapping function without model selection for a model of class 'gamMRSea'
#'
#' This fuction performs a specified number of bootstrapping iterations using CReSS/SALSA for fitting the count model. See below for details.
#'
#' @param predictionGrid The prediction grid data
#' @param model.obj The best fitting \code{CReSS} model for the original count data. Should be geeglm or a Poisson/Binomial GLM (not quasi).
#' @param splineParams The object describing the parameters for fitting the one and two dimensional splines
#' @param g2k (N x k) matrix of distances between all prediction points (N) and all knot points (k)
#' @param rename A vector of column names for which a new column needs to be created for the bootstrapped data. This defaults to \code{segment.id} for line transects (which is required for \code{create.bootcount.data}), others might be added.
#' A new column with new ids will automatically be created for the column listed in \code{resample}. In case of nearshore data, this argument is ignored.
#' @param B Number of bootstrap iterations
#' @param name Analysis name. Required to avoid overwriting previous bootstrap results. This name is added at the beginning of "predictionboot.RData" when saving bootstrap predictions.
#' @param seed Set the seed for the bootstrap sampling process.
#' @param nCores Set the number of computer cores for the bootstrap process to use (default = 1).  The more cores the faster the proces but be wary of over using the cores on your computer. If \code{nCores} > (number of computer cores - 2), the function defaults to \code{nCores} = (number of computer cores - 2).  Note: On a Mac computer the parallel code does not compute so use nCores=1.
#'
#' @details
#'
#' The following steps are performed for each iteration:
#'
#' - coefficients are resampled from a multivariate normal distribution defined by MLE and COV from the best fitting count model
#'
#' - predictions are made to the prediction data using the resampled coefficients
#'
#' @return
#' The function returns a matrix of bootstrap predictions. The number of rows is equal to the number of rows in predictionGrid.  The number of columns is equal to \code{B}.  The matrix may be very large and so is stored directly into the working directory as a workspace object: '"name"predictionboot.RObj'.  The object inside is called \code{bootPreds}.
#'
#'
#'
#' @export

do.bootstrap.cress.robust<-function(model.obj, predictionGrid, splineParams=NULL, g2k=NULL, B, robust=T,name=NULL, seed=12345, nCores=1, cat.message=TRUE){

  require(Matrix)
  require(mvtnorm)

  #require(mvtnorm)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~ PARALLEL CODE ~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(nCores>1){
    #require(parallel)
    if(cat.message) {cat('Code in parallel: No progress guide printed')}

    computerCores <- getOption("cl.cores", detectCores())
    if(nCores>(computerCores-2)){nCores<-(computerCores-2)}
    myCluster <- makeCluster(nCores) ; myCluster # initialise the cluster
    clusterExport(myCluster, ls(), envir=environment()) # export all objects to each cluster
    # export directory and functions to each cluster
    clusterEvalQ(myCluster, {
      #setwd('../Dropbox/SNH/Nhats/')
      #setwd(dir())
      require(mrds)
      require(MRSea)
      require(geepack)
      require(mvtnorm)
      require(mgcv)
      require(splines)
    })

    # only do parametric boostrap if no data re-sampling and no nhats provided
    Routputs<-parLapply(myCluster, 1:B, function(b){

      set.seed(seed + b)
      # setting up CReSS related variables
      attributes(model.obj$formula)$.Environment<-environment()
      radii<-splineParams[[1]]$radii
      d2k<-splineParams[[1]]$dist
      radiusIndices<-splineParams[[1]]$radiusIndices
      aR<- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]

      est<- coefficients(model.obj)
      if(robust==T){
        vbeta<-summary(model.obj)$cov.robust
        samplecoeff<-NULL
        try(samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd')))
        if(is.null(samplecoeff)){
          vbeta<-as.matrix(nearPD(as.matrix(summary(model.obj)$cov.robust))$mat)
          samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd'))
        }
      }
      else{
        vbeta<-summary(model.obj)$cov.scaled
        samplecoeff<-NULL
        try(samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd')))
        if(is.null(samplecoeff)){
          vbeta<-as.matrix(nearPD(as.matrix(summary(model.obj)$cov.scaled))$mat)
          samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd'))
        }
      }
      # make predictions
      dists<- g2k
      rpreds<-predict.gamMRSea(newdata=predictionGrid, g2k=dists, object=model.obj, type='response', coeff=samplecoeff)
      bootPreds<- rpreds
      return(bootPreds)
    })

    stopCluster(myCluster)
    bootPreds<- matrix(unlist(Routputs), ncol=B)

  } # end ncores >1 loop

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~ SINGLE CORE CODE ~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(nCores==1){
    # object for storing the predictions
    bootPreds<- matrix(NA, nrow=nrow(predictionGrid), ncol=B)
    if(cat.message) {cat('Not in parallel: progress guide printed below \n')}

    attributes(model.obj$formula)$.Environment<-environment()
    radii<-splineParams[[1]]$radii
    d2k<-splineParams[[1]]$dist
    radiusIndices<-splineParams[[1]]$radiusIndices
    aR<- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]

    for(b in 1:B){#
      set.seed(seed + b)
      #cat('Parametric boostrp only')

      if(cat.message) {if((b/10)%%1 == 0){cat(b, '\n')}else{cat('.')}}

      # sample from multivariate normal
      est<- coefficients(model.obj)
      if(robust==T){
        vbeta<-summary(model.obj)$cov.robust
        samplecoeff<-NULL
        try(samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd')))
        if(is.null(samplecoeff)){
          vbeta<-as.matrix(nearPD(as.matrix(summary(model.obj)$cov.robust))$mat)
          samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd'))
        }
      }
      else{
        vbeta<-summary(model.obj)$cov.scaled
        samplecoeff<-NULL
        try(samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd')))
        if(is.null(samplecoeff)){
          vbeta<-as.matrix(nearPD(as.matrix(summary(model.obj)$cov.scaled))$mat)
          samplecoeff<- as.numeric(rmvnorm(1,est,vbeta, method='svd'))
        }
      }
      # make predictions
      dists<- g2k
      rpreds<-predict.gamMRSea(newdata=predictionGrid, g2k=dists, object=model.obj, type='response', coeff=samplecoeff)
      bootPreds[,b]<- rpreds
    }

  } # end ncores=1 loop

  # save predictions  -  and other data? e.g. parameter values?
  #save(bootPreds, file=paste(name,"predictionboot.RData",sep=""), compress='bzip2')
return(bootPreds)
} # end of function

