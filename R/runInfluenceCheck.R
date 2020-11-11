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
  require(mgcv)
  inflStore <- matrix(0, nrow = length(unique(idUse)), ncol = (length(coef(model)) + 
                                                                 2))
  timeForOne <- system.time(for (i in unique(idUse)) {
    rowsToDel <- which(id == i)
    pos <- which(i == idUse)
    newData <- dat[-rowsToDel, ]
    
    presPred <- as.matrix(model.matrix(model)[rowsToDel, ]) %*% inflStore[pos, 1:length(coef(model))]
    inflStore[pos, ncol(inflStore)] <- sum((response[rowsToDel] - family(model)$linkinv(presPred))^2)
    model.det <- det(summary(model)$cov.scaled)
    
    if ("splineParams" %in% names(model)) {
      orig.dist<-model$splineParams[[1]]$dist
      model$splineParams[[1]]$dist<- model$splineParams[[1]]$dist[-rowsToDel,]
    }
    
    options(warn = -1)
    newMod <- update(model, . ~ ., data = newData)
    if("panels" %in% names(newMod)){
      newMod$panels<-newMod$panels[-rowsToDel]
    }
    
    inflStore[pos, 1:length(coef(newMod))] <- newMod$coefficients
    inflStore[pos, (ncol(inflStore) - 1)] <- det(summary(newMod)$cov.scaled)/model.det
  })
  print(paste("Calculating the influence measures will take approximately ", 
              round(timeForOne[3] * length(unique(id))/60, 0), " minutes", 
              sep = ""))
  options(warn = 0)
}


