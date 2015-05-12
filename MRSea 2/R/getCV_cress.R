#-----------------------------------------------------------------------------
#' Calculate cross-validation score for a CReSS type model
#'
#' @param data Data frame containing columns of covariates contained in \code{baseModel}.
#' @param baseModel glm or CReSS type model object
#' @param splineParams list object containing information for fitting one and two dimensional splines. See \code{\link{makesplineParams}} for more details.
#' @param vector Logical indicating whether the vector of scores is returned (TRUE) or if the mean score is returned (FALSE).
#' 
#' @details
#' There must be a column in the data called \code{foldid}, which can be created using \code{\link{getCVids}}.  This column defines the folds of data for the CV calculation.
#' 
#' The cost function for this CV is a mean squared error.
#' 
#' @examples
#' # load data
#' data(ns.data.re)
#' # load prediction data
#' data(ns.predict.data.re)
#' 
#' splineParams<-makesplineParams(data=ns.data.re, varlist=c('observationhour', 'DayOfMonth'))
#' # set some input info for SALSA
#' ns.data.re$response<- ns.data.re$birds
#' salsa1dlist<-list(fitnessMeasure = 'QICb', minKnots_1d=c(2,2), maxKnots_1d = c(20, 20), 
#'                startKnots_1d = c(2,2), degree=c(2,2), maxIterations = 10, gaps=c(1,1))
#'      
#' # set initial model without the spline terms in there 
#' # (so all other non-spline terms)
#' initialModel<- glm(response ~ as.factor(floodebb) + as.factor(impact) + offset(log(area)), 
#'                   family='quasipoisson',data=ns.data.re)
#'                   
#' # run SALSA
#' salsa1dOutput<-runSALSA1D(initialModel, salsa1dlist, varlist=c('observationhour','DayOfMonth'), 
#'             factorlist=c('floodebb', 'impact'), ns.predict.data.re, splineParams=splineParams)
#'  
#' # make blocking structure and fold structure
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
#'                     ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)
#' ns.data.re$foldid<-getCVids(ns.data.re, folds=5, block='blockid')
#' 
#' # calculate CV
#' cv1<-getCV_CReSS(ns.data.re, salsa1dOutput$bestModel, salsa1dOutput$splineParams)
#' 
#' @export
#' 
getCV_CReSS<-function(datain, baseModel, splineParams, vector=FALSE){

  #data<-baseModel$data
  
#   if(is.null(block)==FALSE){
#     data$foldid<-data[,block]  
#   }
  
  d2k<-splineParams[[1]]$dist
  radiusIndices <-splineParams[[1]]$radiusIndices
  radii <- splineParams[[1]]$radii
  aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
  
  attributes(baseModel$formula)$.Environment<-environment()
  # calculate cross-validation
  nfolds<-length(unique(datain$foldid))
  dists<-d2k
  store<- matrix(0, nrow=nfolds, ncol=1)
  for(f in 1:nfolds){
    dists<-d2k[datain$foldid!=f,]
    
    eval(parse(text=paste(substitute(datain), "=datain[datain$foldid!=", f, ",]", sep="")))
    
    eval(parse(text=paste("foldedFit<- update(baseModel, .~., data=", substitute(datain), ")", sep="")))
    
    if(length(coef(foldedFit))==1){
      dists<-d2k[datain$foldid==f,]
      
      if(is.null(baseModel$offset)){
        predscv<- baseModel$family$linkinv(as.matrix(model.matrix(baseModel)[datain$foldid==f])%*%coef(foldedFit))
      }else{
        predscv<- baseModel$family$linkinv(as.matrix(model.matrix(baseModel)[datain$foldid==f])%*%coef(foldedFit)) * baseModel$family$linkinv(baseModel$offset)[datain$foldid==f]  
      }
      
    }else{  
      dists<-d2k[datain$foldid==f,]
      
      if(is.null(baseModel$offset)){
      predscv<- baseModel$family$linkinv(model.matrix(baseModel)[datain$foldid==f,]%*%coef(foldedFit))
      }else{
        predscv<- baseModel$family$linkinv(model.matrix(baseModel)[datain$foldid==f,]%*%coef(foldedFit)) * baseModel$family$linkinv(baseModel$offset)[datain$foldid==f]
      }
      
    }        
    store[f]<- mean((datain$response[datain$foldid==f]-predscv)**2)
  }
  if(vector==TRUE){
    return(store)
  }else{
    return(mean(store))
  }
}

