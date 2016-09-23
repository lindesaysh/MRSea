#' Function to make model of class \code{gamMRSea}
#'
#'
#' @param model
#' @param panelid
#' @param splineParams
#'
#'
#' @example
#'

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
