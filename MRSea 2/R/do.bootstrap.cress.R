#' Bootstrapping function without model selection using CReSS/SALSA for fitting the second stage count model
#' 
#' This fuction performs a specified number of bootstrapping iterations using CReSS/SALSA for fitting the 
#' second stage count model. See below for details. 
#' 
#' @param orig.data The original data. In case \code{ddf.obj} is specified, this should be the original distance data. In case \code{ddf.obj} is \code{NULL}, it should have the format equivalent to \code{count.data} where each record represents the summed up counts at the segments. 
#' @param predict.data The prediction grid data 
#' @param ddf.obj The ddf object created for the best fitting detection model. Defaults to \code{NULL} for when nodetection function object available. 
#' @param model.obj The best fitting \code{CReSS} model for the original count data
#' @param splineParams The object describing the parameters for fitting the one and two dimensional splines
#' @param g2k (N x k) matrix of distances between all prediction points (N) and all knot points (k)
#' @param resample Specifies the resampling unit for bootstrapping, default is \code{transect.id}. Must match a column name in \code{dis.data} exactly
#' @param rename A vector of column names for which a new column needs to be created for the bootstrapped data. This defaults to \code{segment.id} for line transects (which is required for \code{create.bootcount.data}), others might be added. 
#' A new column with new ids will automatically be created for the column listed in \code{resample}. In case of nearshore data, this argument is ignored. 
#' @param stratum The column name in \code{orig.data} that identifies the different strata. The default \code{NULL} returns un-stratified bootstrap data. In case of nearshore data, this argument is ignored.  
#' @param B Number of bootstrap iterations
#' @param name Analysis name. Required to avoid overwriting previous bootstrap results. This name is added at the beginning of "predictionboot.RData" when saving bootstrap predictions. 
#' @param save.data If TRUE, all created bootstrap data will be saved as an RData object in the working directory at each iteration, defaults to FALSE
#' @param nhats (default = FALSE).  If you have calculated bootstrap NHATS because there is no simple ddf object then a matrix of these may be fed into the function.  The number of columns of data should >= B.  The rows must be equal to those in \code{orig.data} and \code{d2k} and \emph{must} be in matching order.
#' @param seed Set the seed for the bootstrap sampling process.
#' @param nCores Set the number of computer cores for the bootstrap process to use (default = 1).  The more cores the faster the proces but be wary of over using the cores on your computer. If \code{nCores} > (number of computer cores - 2), the function defaults to \code{nCores} = (number of computer cores - 2).  Note: On a Mac computer the parallel code does not compute so use nCores=1.
#' 
#' @details
#' In case of distance sampling data, the following steps are performed for each iteration: 
#' 
#' - the original data is bootstrapped
#' 
#' - a detection function is fitted to the bootstrapped data
#' 
#' - a count model is fitted to the bootstrapped data
#' 
#' - coefficients are resampled from a multivariate normal distribution defined by MLE and COV from count model
#' 
#' - predictions are made to the prediction data using the resampled coefficients 
#' 
#' In case of count data, the following steps are performed for each iteration:
#' 
#' - coefficients are resampled from a multivariate normal distribution defined by MLE and COV from the best fitting count model
#' 
#' - predictions are made to the prediction data using the resampled coefficients 
#' 
#' @return
#' The function returns a matrix of bootstrap predictions. The number of rows is equal to the number of rows in predict.data.  The number of columns is equal to \code{B}.  The matrix may be very large and so is stored directly into the working directory as a workspace object: '"name"predictionboot.RObj'.  The object inside is called \code{bootPreds}.
#' 
#' @examples
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # offshore redistribution data
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' data(dis.data.re)
#' data(predict.data.re)
#' data(knotgrid.off)

#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # distance sampling
#' dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
#' result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
#'         meta.data=list(width=250))
#' dis.data.re<-create.NHAT(dis.data.re,result)
#' count.data<-create.count.data(dis.data.re)
#' 
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # spatial modelling
#' splineParams<-makesplineParams(data=count.data, varlist=c('depth'))
#' #set some input info for SALSA
#' count.data$response<- count.data$NHAT
#' # make distance matrices for datatoknots and knottoknots
#' distMats<-makeDists(cbind(count.data$x.pos, count.data$y.pos), na.omit(knotgrid.off))
#' # choose sequence of radii
#' r_seq<-getRadiiChoices(8,distMats$dataDist)
#' # set initial model without the spatial term
#' initialModel<- glm(response ~ as.factor(season) + as.factor(impact) + offset(log(area)),  
#'                 family='quasipoisson', data=count.data)
#' # make parameter set for running salsa2d
#' salsa2dlist<-list(fitnessMeasure = 'QICb', knotgrid = knotgrid.off, knotdim=c(26,14), startKnots=4, minKnots=4, 
#'                  maxKnots=20, r_seq=r_seq, gap=4000, interactionTerm="as.factor(impact)")
#' salsa2dOutput_k6<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist, 
#'                    k2k=distMats$knotDist, splineParams=splineParams) 
#' 
#' splineParams<-salsa2dOutput_k6$splineParams
#' # specify parameters for local radial function:
#' radiusIndices <- splineParams[[1]]$radiusIndices
#' dists <- splineParams[[1]]$dist
#' radii <- splineParams[[1]]$radii
#' aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
#' count.data$blockid<-paste(count.data$transect.id, count.data$season, count.data$impact, sep='')
#' # Re-fit the chosen model as a GEE (based on SALSA knot placement) and GEE p-values
#' geeModel<- geeglm(formula(salsa2dOutput_k6$bestModel), data=count.data, family=poisson, id=blockid)
#' dists<-makeDists(cbind(predict.data.re$x.pos, predict.data.re$y.pos), na.omit(knotgrid.off), 
#'        knotmat=FALSE)$dataDist
#' 
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # bootstrapping
#' do.bootstrap.cress(dis.data.re, predict.data.re, ddf.obj=result, geeModel, splineParams, 
#'               g2k=dists, resample='transect.id', rename='segment.id', stratum='survey.id', 
#'               B=4, name="cress", save.data=FALSE, nhats=FALSE, nCores=1)
#' load("cresspredictionboot.RData") # loading the bootstrap predictions into the workspace
#' # look at the first 6 lines of the bootstrap predictions (on the scale of the response)
#' head(bootPreds) 
#' 
#' \dontrun{
#' # In parallel (Note: windows machines only)
#' require(parallel)
#' do.bootstrap.cress(dis.data.re, predict.data.re, ddf.obj=result, geeModel, splineParams,
#'                 g2k=dists, resample='transect.id', rename='segment.id', stratum='survey.id',
#'                 B=4, name="cress", save.data=FALSE, nhats=FALSE, nCores=4)
#' load("cresspredictionboot.RData") # loading the bootstrap predictions into the workspace
#' # look at the first 6 lines of the bootstrap predictions (on the scale of the response)
#' head(bootPreds)
#' }
#' 
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # nearshore redistribution data
#' # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' \dontrun{
#' do.bootstrap.cress(ns.data.re, ns.predict.data.re, ddf.obj=NULL, geeModel, splineParams, 
#'              g2k=dists, resample='transect.id', rename='segment.id', stratum=NULL, 
#'              B=2, name="cress", save.data=FALSE, nhats=FALSE)
#' load("cresspredictionboot.RData") # loading the predictions into the workspace
#' # look at the first 6 lines of the bootstrap predictions (on the scale of the response)
#' head(bootPreds)}
#' 
#' @export

do.bootstrap.cress<-function(orig.data,predict.data,ddf.obj=NULL,model.obj,splineParams, g2k, resample="transect.id",rename="segment.id",stratum=NULL,B,name=NULL, save.data=FALSE, nhats=FALSE, seed=12345, nCores=1){
  
  if(is.null(nhats)){nhats<-FALSE}
  #require(mvtnorm)
  if(is.null(ddf.obj)==F){
    # fixing the settings for the detection model
    key <- ddf.obj$ds$aux$ddfobj$type
    scale.formula<-ddf.obj$ds$aux$ddfobj$scale$formula
    point<-ddf.obj$ds$aux$point
    binned<-ddf.obj$ds$aux$binned
    breaks<-ddf.obj$ds$aux$breaks
    width<-ddf.obj$ds$aux$width
    my.formula=eval(parse(file="",text=paste("~mcds(key = '",key,"', formula=",scale.formula,")",sep="")))
    #result=ddf(dsmodel = as.formula(my.formula), data = data$dis.object, method = "ds", meta.data = list(point = point, width = data$w, binned = data$binned))
  }
  
  
  if(is.null(ddf.obj)==F){print('Bootstrap with data resampling')
  }else{print('Bootstrap without data resampling')}
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~ PARALLEL CODE ~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(nCores>1){
    #require(parallel)
    cat('Code in parallel: No progress guide printed')
    
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
    })
    
    # only do parametric boostrap if no data re-sampling and no nhats provided
    if(is.null(ddf.obj) & nhats==FALSE){
      
      
      Routputs<-parLapply(myCluster, 1:B, function(b){
        
        set.seed(seed + b)
        # setting up CReSS related variables
        attributes(model.obj$formula)$.Environment<-environment()
        radii<-splineParams[[1]]$radii
        d2k<-splineParams[[1]]$dist
        radiusIndices<-splineParams[[1]]$radiusIndices
        aR<- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
        
        # sample from multivariate normal
        est<- coefficients(model.obj)
        varcov<-as.matrix(summary(model.obj)$cov.unscaled)
        samplecoeff<-NULL
        try(samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd')))
        if(is.null(samplecoeff)){
          varcov<-as.matrix(nearPD(vcov(model.obj))$mat)
          samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd'))
        }
        
        # make predictions
        dists<- g2k
        rpreds<-predict.cress(predict.data, splineParams, dists, model.obj, type='response', coeff=samplecoeff)
        bootPreds<- rpreds
        return(bootPreds)
      })
      
      stopCluster(myCluster)
      bootPreds<- matrix(unlist(Routputs), ncol=B)
      
    }else{
      Routputs<-parLapply(myCluster, 1:B, function(b){
        
        set.seed(seed + b)
      
        # setting up CReSS related variables
        attributes(model.obj$formula)$.Environment<-environment()
        radii<-splineParams[[1]]$radii
        d2k<-splineParams[[1]]$dist
        radiusIndices<-splineParams[[1]]$radiusIndices
        aR<- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
        
        if(is.null(ddf.obj)==F){
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~ DISTANCE BOOTSTRAP  ~~~~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          # create bootstrap data
          bootstrap.data<-create.bootstrap.data(orig.data, stratum = stratum)
          
          if (save.data==T){
            filename=paste("bootstrap.data_",b,".RData",sep="")
            save(bootstrap.data,file=filename)
          }
          
          # fit the same detection model as in ddf.obj to the bootstrap data
          result.boot<-ddf(dsmodel=as.formula(my.formula), data=bootstrap.data, method="ds", meta.data=list(width=width, binned=binned, breaks=breaks))
          
          # create count data
          bootstrap.data<-create.NHAT(bootstrap.data, result.boot)
          boot.count.data<-create.bootcount.data(bootstrap.data)
          #if(is.null(boot.count.data$response)){stop('No column in orig.data called response')}
          boot.count.data$response<-boot.count.data$NHAT
          boot.count.data$transect.id2<-as.factor(boot.count.data$transect.id2)
          
          # sort by segment.id2 and then update dist accordingly
          gamsort<-boot.count.data[order(boot.count.data$segment.id),]
          newid<-NULL
          for(n in 1:length(unique(boot.count.data$transect.id2))){
            newid<-c(newid, which(match(x=gamsort$transect.id2, table=unique(gamsort$transect.id2)[n])==1))
          }
          gamsort<-gamsort[newid,]
          boot.count.data<-gamsort
          
          # ******** Fit models ********** #
          #print(paste('zeroblocks:', length(which(tapply(boot.count.data$response, boot.count.data$transect.id2, sum)==0)), sep=''))
          
          dists<- d2k[boot.count.data$segment.id,]
          #geefit<- geeglm(model.obj$formula, id=transect.id2, family=model.obj$family$family, data=boot.count.data)
          geefit<- eval(parse(text=paste("geeglm(model.obj$formula, id=transect.id2, family=", model.obj$family$family,"(link='", model.obj$family$link, "'), data=boot.count.data)", sep='')))
          
          # sample from multivariate normal
          est<- coefficients(geefit)
          varcov<-as.matrix(summary(geefit)$cov.unscaled)
          samplecoeff<-NULL
          try(samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd')), silent=T)
          if(is.null(samplecoeff)){
            varcov<-as.matrix(nearPD(vcov(model.obj))$mat)
            samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd'))
            #cat('*')
          }
          
          # make predictions
          dists<- g2k
          rpreds<-predict.cress(predict.data, splineParams, dists, geefit, type='response', coeff=samplecoeff)
          #x2<- data.frame(response=rep(1,nrow(predict.data)), predict.data)
          #fakemodel<- glm(geefit$formula, data=x2, family=quasipoisson)
          #rpreds<- exp(model.matrix(fakemodel)%*%samplecoeff)* predict.data$area
          bootPreds<- rpreds
          
        }else{
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~ NO DISTANCE BOOTSTRAP  ~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          #orig.data$response<-orig.data$count
          if(is.null(orig.data$response)){stop('No column in orig.data called response')}
          
          orig.data$response<-nhats[,b]
        
          # ******** Fit models ********** #
          
          dists<- d2k
          geefit<- eval(parse(text=paste("geeglm(model.obj$formula, id=blockid, family=", model.obj$family$family,"(link='", model.obj$family$link, "'), data=orig.data)", sep='')))
          
          # sample from multivariate normal
          est<- coefficients(geefit)
          varcov<-as.matrix(summary(geefit)$cov.unscaled)
          samplecoeff<-NULL
          try(samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd')))
          if(is.null(samplecoeff)){
            varcov<-as.matrix(nearPD(vcov(geefit))$mat)
            samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd'))
          }
          
          # make predictions
          dists<- g2k
          rpreds<-predict.cress(predict.data, splineParams, dists, geefit, type='response', coeff=samplecoeff)
          bootPreds<- rpreds
        }  # end of else statement
        
        return(bootPreds)
        # end of cluster for loop
      })
      
      stopCluster(myCluster)
      
      bootPreds<- matrix(unlist(Routputs), ncol=B)
    } # end else statement
  
  } # end ncores >1 loop
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~ SINGLE CORE CODE ~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(nCores==1){
    # object for storing the predictions
    bootPreds<- matrix(NA, nrow=nrow(predict.data), ncol=B)
    cat('Not in parallel: progress guide printed below \n')
    
    attributes(model.obj$formula)$.Environment<-environment()
    radii<-splineParams[[1]]$radii
    d2k<-splineParams[[1]]$dist
    radiusIndices<-splineParams[[1]]$radiusIndices
    aR<- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
    
    
    if(is.null(ddf.obj) & nhats==FALSE){
      
      for(b in 1:B){#
        set.seed(seed + b)
        #cat('Parametric boostrp only')
        
        if((b/10)%%1 == 0){cat(b, '\n')}else{cat('.')}
        
        # sample from multivariate normal
        est<- coefficients(model.obj)
        varcov<-as.matrix(summary(model.obj)$cov.unscaled)
        samplecoeff<-NULL
        try(samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd')))
        if(is.null(samplecoeff)){
          varcov<-as.matrix(nearPD(vcov(model.obj))$mat)
          samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd'))
        }
        # make predictions
        dists<- g2k
        rpreds<-predict.cress(predict.data, splineParams, dists, model.obj, type='response', coeff=samplecoeff)
        bootPreds[,b]<- rpreds
      }
    }else{
      for (b in 1:B){
        
        set.seed(seed + b)
        
        if((b/10)%%1 == 0){cat(b, '\n')}else{cat('.')}
        
        if(is.null(ddf.obj)==F){
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~ DISTANCE BOOTSTRAP  ~~~~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # create bootstrap data
          bootstrap.data<-create.bootstrap.data(orig.data, stratum = stratum)
          
          if (save.data==T){
            filename=paste("bootstrap.data_",b,".RData",sep="")
            save(bootstrap.data,file=filename)
          }
          
          # fit the same detection model as in ddf.obj to the bootstrap data
          result.boot<-ddf(dsmodel=as.formula(my.formula), data=bootstrap.data, method="ds", meta.data=list(width=width, binned=binned, breaks=breaks))
          
          # create count data
          bootstrap.data<-create.NHAT(bootstrap.data, result.boot)
          boot.count.data<-create.bootcount.data(bootstrap.data)
          #if(is.null(boot.count.data$response)){stop('No column in orig.data called response')}
          boot.count.data$response<-boot.count.data$NHAT
          boot.count.data$transect.id2<-as.factor(boot.count.data$transect.id2)
          
          # sort by segment.id2 and then update dist accordingly
          gamsort<-boot.count.data[order(boot.count.data$segment.id),]
          newid<-NULL
          for(n in 1:length(unique(boot.count.data$transect.id2))){
            newid<-c(newid, which(match(x=gamsort$transect.id2, table=unique(gamsort$transect.id2)[n])==1))
          }
          gamsort<-gamsort[newid,]
          boot.count.data<-gamsort
          
          # ******** Fit models ********** #
          #print(paste('zeroblocks:', length(which(tapply(boot.count.data$response, boot.count.data$transect.id2, sum)==0)), sep=''))
          
          dists<- d2k[boot.count.data$segment.id,]
          #geefit<- geeglm(model.obj$formula, id=transect.id2, family=model.obj$family$family, data=boot.count.data)
          geefit<- eval(parse(text=paste("geeglm(model.obj$formula, id=transect.id2, family=", model.obj$family$family,"(link='", model.obj$family$link, "'), data=boot.count.data)", sep='')))
          
          # sample from multivariate normal
          est<- coefficients(geefit)
          varcov<-as.matrix(summary(geefit)$cov.unscaled)
          samplecoeff<-NULL
          try(samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd')), silent=T)
          if(is.null(samplecoeff)){
            varcov<-as.matrix(nearPD(vcov(model.obj))$mat)
            samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd'))
            cat('*')
          }
          
          # make predictions
          dists<- g2k
          rpreds<-predict.cress(predict.data, splineParams, dists, geefit, type='response', coeff=samplecoeff)
          #x2<- data.frame(response=rep(1,nrow(predict.data)), predict.data)
          #fakemodel<- glm(geefit$formula, data=x2, family=quasipoisson)
          #rpreds<- exp(model.matrix(fakemodel)%*%samplecoeff)* predict.data$area
          bootPreds[,b]<- rpreds
          
        }else{
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~ NO DISTANCE BOOTSTRAP  ~~~~~~~~~~~~~~~~
          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          #orig.data$response<-orig.data$count
          if(is.null(orig.data$response)){stop('No column in orig.data called response')}
          
          #if(nhats==TRUE){
          orig.data$response<-nhats[,b]
          #}
          
          # ******** Fit models ********** #
          
          dists<- d2k
          #geefit<- geeglm(model.obj$formula, id=blockid, family=model.obj$family$family, data=orig.data)
          geefit<- eval(parse(text=paste("geeglm(model.obj$formula, id=blockid, family=", model.obj$family$family,"(link='", model.obj$family$link, "'), data=orig.data)", sep='')))
          
          # sample from multivariate normal
          est<- coefficients(geefit)
          varcov<-as.matrix(summary(geefit)$cov.unscaled)
          samplecoeff<-NULL
          try(samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd')))
          if(is.null(samplecoeff)){
            varcov<-as.matrix(nearPD(vcov(geefit))$mat)
            samplecoeff<- as.numeric(rmvnorm(1,est,varcov, method='svd'))
          }
          
          # make predictions
          dists<- g2k
          rpreds<-predict.cress(predict.data, splineParams, dists, geefit, type='response', coeff=samplecoeff)
          #x2<- data.frame(response=rep(1,nrow(predict.data)), predict.data)
          #fakemodel<- glm(geefit$formula, data=x2, family=quasipoisson)
          #rpreds<- exp(model.matrix(fakemodel)%*%samplecoeff)* predict.data$area
          bootPreds[,b]<- rpreds
        }  # end of else statement
        
      } # end of for loop
    } # end of else statement
  } # end ncores=1 loop
  
  # save predictions  -  and other data? e.g. parameter values?
  save(bootPreds, file=paste(name,"predictionboot.RData",sep=""), compress='bzip2')
  
} # end of function

