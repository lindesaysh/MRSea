#-----------------------------------------------------------------------------
#' Calculate cross-validation score for a CReSS type model
#'
#' @param datain Data frame containing columns of covariates contained in \code{baseModel}.
#' @param baseModel 'glm' or 'gamMRSea' model object
#' @param splineParams list object containing information for fitting one and two dimensional splines. See \code{\link{makesplineParams}} for more details.
#' @param vector Logical indicating whether the vector of scores is returned (TRUE) or if the mean score is returned (FALSE).
#' 
#' @details
#' There must be a column in the data called \code{foldid}, which can be created using \code{\link{getCVids}}.  This column defines the folds of data for the CV calculation.  If this is not provided, by default random allocation is given to 5 folds.
#' 
#' The cost function for this CV is a mean squared error.
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' ns.data.re$foldid<-getCVids(ns.data.re, folds=5)
#'  
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#'               
#' # calculate CV
#' getCV_CReSS(ns.data.re, model)
#' 
#' @author LAS Scott-Hayward (University of St Andrews)
#' @export
#' 
getCV_CReSS<-function(datain, baseModel, splineParams=NULL, vector=FALSE){

  #data<-baseModel$data
  
#   if(is.null(block)==FALSE){
#     data$foldid<-data[,block]  
#   }
  
  if(is.null(datain$foldid)){
    datahere<-datain
    datahere$foldid<-getCVids(datain, folds = 5)
    warning("Fold ID not specified as column in the data.  Random allocation to 5 folds used here by default.  To change, use getCVids() to create your own foldid column in the data set.")
  }else{
    datahere<-datain
  }
  
  if(is.null(splineParams)){
    splineParams<-baseModel$splineParams  
  }
  # 
   d2k<-splineParams[[1]]$dist
  # radiusIndices <-splineParams[[1]]$radiusIndices
  # radii <- splineParams[[1]]$radii
  # aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]

  
  
  attributes(baseModel$formula)$.Environment<-environment()
  # calculate cross-validation
  nfolds<-length(unique(datahere$foldid))
  #dists<-d2k
  store<- matrix(0, nrow=nfolds, ncol=1)
  for(f in 1:nfolds){
    baseModel$splineParams[[1]]$dist<-d2k[datahere$foldid!=f,]
    splineParams[[1]]$dist<-d2k[datahere$foldid!=f,]
    
    eval(parse(text=paste(substitute(datain), "=datahere[datahere$foldid!=", f, ",]", sep="")))
    
    eval(parse(text=paste("foldedFit<- update(baseModel, .~., data=", substitute(datain), ")", sep="")))
    # one parameter model
    if(length(coef(foldedFit))==1){
      #dists<-d2k[datahere$foldid==f,]
      splineParams[[1]]$dist<-d2k
      baseModel$splineParams[[1]]$dist<-d2k
      
      #offset or not?
      if(is.null(baseModel$offset)){
        predscv<- baseModel$family$linkinv(as.matrix(model.matrix(baseModel)[datahere$foldid==f])%*%coef(foldedFit))
      }else{
        predscv<- baseModel$family$linkinv(as.matrix(model.matrix(baseModel)[datahere$foldid==f])%*%coef(foldedFit)) * baseModel$family$linkinv(baseModel$offset)[datahere$foldid==f]  
      }
      
    }else{  
      #dists<-d2k[datahere$foldid==f,]
      splineParams[[1]]$dist<-d2k
      baseModel$splineParams[[1]]$dist<-d2k
      
      if(is.null(baseModel$offset)){
      predscv<- baseModel$family$linkinv(model.matrix(baseModel)[datahere$foldid==f,]%*%coef(foldedFit))
      }else{
        predscv<- baseModel$family$linkinv(model.matrix(baseModel)[datahere$foldid==f,]%*%coef(foldedFit)) * baseModel$family$linkinv(baseModel$offset)[datahere$foldid==f]
      }
      
    }  
    # check for response variable
    if(nrow(baseModel$data)!=length(baseModel$model[[1]])){
      props<-datahere$successes[datahere$foldid==f]/(datahere$successes[datahere$foldid==f] + datahere$failures[datahere$foldid==f])
      store[f]<- mean((props-predscv)**2)
    }else{
      store[f]<- mean((baseModel$y[datahere$foldid==f]-predscv)**2)
    }
    
  }
  if(vector==TRUE){
    return(store)
  }else{
    return(mean(store))
  }
}

