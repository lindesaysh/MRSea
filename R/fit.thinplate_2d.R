#' Function to fit a local radial basis function (CReSS) as a two dimensional smooth
#'
#'
#' @export
#'


"fit.thinPlate_2d" <- function(fitnessMeasure, dists,aR,radii,baseModel,radiusIndices,models, currentFit, interactionTerm, data, initDisp, cv.opts, basis='gaussian') {

  if (isS4(baseModel)){
    attributes(baseModel@misc$formula)$.Environment<-environment()
    baseModel@splineParams[[1]]$knotPos<-aR
    baseModel@splineParams[[1]]$radiusIndices<-radiusIndices
  } else {
    attributes(baseModel$formula)$.Environment<-environment()
    baseModel$splineParams[[1]]$knotPos<-aR
    #baseModel$splineParams[[1]]$knotPos<-baseModel$splineParams[[1]]$mapInd[aR]
    baseModel$splineParams[[1]]$radiusIndices<-radiusIndices
  }

  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Fitting Model...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #aR<- aR
  #radiusIndices<- radiusIndices

  #cat("R Indices: ", radiusIndices, "\n")


  # print(aR)
  if(basis=='gaussian'){
    bspl<- "LRF.g(radiusIndices, dists, radii, aR)"  
  }
  
  if(basis=='exponential'){
    bspl<- "LRF.e(radiusIndices, dists, radii, aR)"  
  }

  if(is.null(interactionTerm)){
    test<-paste("update(baseModel, .  ~ . + ",bspl, ")",sep="")
    currentModel<-eval(parse(text=test))
  }else{
    test<-paste("update(baseModel, .  ~  . + ",bspl, "*",interactionTerm, ")", sep="")
    currentModel<-eval(parse(text=test))
  }

  tempFit <- get.measure_2d(fitnessMeasure, currentFit, currentModel,data, dists,aR,radii,radiusIndices, initDisp, cv.opts)$fitStat
  # if(tempFit <= (currentFit+10)){
  #   models[[length(models)+1]] = list(aR,radiusIndices, radii, tempFit)
  # }
  models<-NULL
  #modelinprogress<<-currentModel

  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Model fitted...")
  #print(paste('disp= ', summary(currentModel)$dispersion, ', num knots: ', length(aR), ', fitstat: ',tempFit,sep=''))
  #print("ooooooooooooooooooooooooooooooooooooooo")
  return(list(currentModel=currentModel,models=models, fitStat=tempFit))

}
