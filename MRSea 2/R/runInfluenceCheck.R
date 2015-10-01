
#-----------------------------------------------------------------------------
#' Timing check to see how long it will take to run \code{runInfluence}.
#' 
#' @param model Fitted model object (glm or gam)
#' @param id blocking structure
#' @param d2k (\code{default=NULL}). (n x k) Matrix of distances between all data points in \code{model} and all valid knot locations.
#' @param splineParams (\code{default=NULL}). List object containng output from runSALSA (e.g. knot locations for continuous covariates). See \code{\link{makesplineParams}} for more details of this object. 
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
#'                      ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)

#' model<-geeglm(birds ~ observationhour + as.factor(floodebb) + as.factor(impact), 
#'                family='poisson', data=ns.data.re, id=blockid)
#' 
#' timeInfluenceCheck(model, ns.data.re$blockid)
#'
#' @export
#' 
timeInfluenceCheck<-function(model, id, d2k=NULL, splineParams=NULL){
  

  attributes(model$formula)$.Environment<-environment()
  response<-model$y
    
  if(class(model)[1]=='geeglm'){
    dat<- model$data
  }
  
  if(class(model)[1]=='gam'){
    dat<-model$model 
  }
  
  idUse<- sample(unique(id),1)
  if(length(which(search()=='package:mgcv'))>0){
    detach("package:mgcv")
  }
  
  require(mgcv)
  
  inflStore<- matrix(0,nrow=length(unique(idUse)), ncol=(length(coef(model))+2))
  
  timeForOne<- system.time(for(i in unique(idUse)){
    rowsToDel<- which(idUse==i)
    pos<- which(i==idUse)
    newData<- dat[-rowsToDel,]
    if(is.null(d2k)==F){dists<- d2k[-rowsToDel,]}
    options(warn=-1)
    newMod<-update(model, .~. ,data=newData)
    inflStore[pos,1:length(coef(model))]<-newMod$coefficients
    presPred<- as.matrix(model.matrix(model)[rowsToDel,])%*%inflStore[pos,1:length(coef(model))]
    inflStore[pos,ncol(inflStore)]<-sum((response[rowsToDel]-family(model)$linkinv(presPred))**2)
    inflStore[pos,(ncol(inflStore)-1)]<-det(summary(newMod)$cov.scaled)/det(summary(model)$cov.scaled)
  })
  print(paste("Calculating the influence measures will take approximately ", round(timeForOne[3]*length(unique(id))/60,0), " minutes", sep=""))
  
  options(warn=0)
  
}

