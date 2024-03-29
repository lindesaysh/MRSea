#' Timing check to see how long it will take to run \code{runInfluence}.
#' 
#' @param model Fitted model object (glm, gamMRSea or gam)
#' @param id blocking structure
#' @examples 
#' # load data
#' data(ns.data.re)
#' 
#' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
#'                      ns.data.re$DayOfMonth, sep='')
#' ns.data.re$blockid<-as.factor(ns.data.re$blockid)

#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#' 
#' timeInfluenceCheck(model, ns.data.re$blockid)
#'
#' @export
#' 
timeInfluenceCheck<-function (model, id) {
  attributes(model$formula)$.Environment <- environment()
  response <- model$y
  if (class(model)[1] == "geeglm" | class(model)[1] == "gamMRSea" | 
      class(model)[1] == "glm") {
    dat <- model$data
  }
  if (class(model)[1] == "gam") {
    dat <- model$model
  }
  idUse <- sample(unique(id), 1)
  if (length(which(search() == "package:mgcv")) > 0) {
    detach("package:mgcv")
  }
  suppressMessages({
    library(mgcv, quietly = TRUE)
  })
  
  inflStore <- matrix(0, nrow = length(unique(id)), ncol = (length(coef(model)) + 2))
  timeForOne <- system.time(for (i in unique(idUse)) {
    #find rows to delete
    rowsToDel<- which(id==i)
    pos<- which(i==unique(id))
    newData<- dat[-rowsToDel,]
    
    # make assessment on original model
    nb<- as.matrix(model.matrix(model)[c(1),])
    nc<- as.matrix(inflStore[pos,1:length(coef(model))])
    
    if(length(pos==1)){
      presPred <- t(nb)%*%nc
    }else{
      presPred <- nb%*%nc  
    }
    
    inflStore[pos,ncol(inflStore)]<-sum((response[rowsToDel]-c(family(model)$linkinv(presPred)))**2)
    
    if(class(model)[1]=='gamMRSea'){
      model.det<-det(summary(model)$cov.robust)
    }else{
      model.det<-det(summary(model)$cov.scaled)
    }
    
    # update model for reduced data size (includes data and distance matrix)
    newMod<-model
    if ("splineParams" %in% names(model)) {
      newMod$splineParams[[1]]$dist<- newMod$splineParams[[1]]$dist[-rowsToDel,]
      splineParams = newMod$splineParams
    }
    
    if(class(model)[1]=='gamMRSea'){
      newpanel<-id[-rowsToDel]
      if(is.factor(newpanel)){
        newpanel<-droplevels(newpanel)
      }
      newMod<-update(newMod, .~. ,data=newData, panels=newpanel, splineParams = splineParams)
      newmod.det<-det(summary(newMod)$cov.robust)
    }else{
      newMod<-update(newMod, .~. ,data=newData)
      newmod.det<-det(summary(newMod)$cov.scaled)
    }
    
    # make assessment on new model
    inflStore[pos,(1:length(coef(newMod)))]<-newMod$coefficients
    inflStore[pos,(ncol(inflStore)-1)]<-newmod.det/model.det
  })
  
  print(paste("Calculating the influence measures will take approximately ", 
              round((timeForOne[3]) * length(unique(id))/60), " minutes", 
              sep = ""))
  
}


