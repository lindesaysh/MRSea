#' Function to make model of class \code{gamMRSea}
#'
#'
#' @param model  model object of class glm or gamMRSea
#' @param panelid vector of length of the data containing the panel identification for each row of data
#' @param splineParams MRSea based list object
#' @param varshortnames vector containing the short names for each variable.  These are used in summary and anova
#' @param gamMRSea logical stating whether the call of the model should be changed to 'gamMRSea'
#'
#'
#' @export

make.gamMRSea<-function(model, panelid=NULL, splineParams=NULL, varshortnames=NULL, gamMRSea=FALSE){
  newmodel<-model
  if(class(model)[1]!='gamMRSea'){
    class(newmodel)<-c('gamMRSea', class(model))
  }
if(is.null(panelid) & is.null(model$panels)){
  newmodel$panels<-1:nrow(model$data)
}
if(!is.null(panelid)){
  newmodel$panels<-panelid
}
if(is.null(panelid) & !is.null(model$panels)){
  newmodel$panels<-model$panels
}

  newmodel$splineParams<-splineParams
  newmodel$varshortnames<-varshortnames

  if(gamMRSea){
    newmodel$call[[1]]<-quote(gamMRSea)
    newmodel$call$splineParams<-quote(splineParams)
  }

  return(newmodel)
}
