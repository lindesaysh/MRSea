#' Function to fit a local radial basis function (CReSS) as a two dimensional smooth
#' 
#' 
#' @export
#' 


"fit.thinPlate_2d" <- function(fitnessMeasure, dists,aR,radii,baseModel,radiusIndices,models, currentFit, interactionTerm, data) {
  
  attributes(baseModel$formula)$.Environment<-environment()
  
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Fitting Model...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #aR<- aR
  #radiusIndices<- radiusIndices
  
  #cat("R Indices: ", radiusIndices, "\n")                         
  
  
  # print(aR)
  bspl<- "LocalRadialFunction(radiusIndices, dists, radii, aR)"
  
  
  if(is.null(interactionTerm)){  
    test<-paste("update(baseModel, .  ~ . + ",bspl, ")",sep="")
    currentModel<-eval(parse(text=test))
  }else{
    test<-paste("update(baseModel, .  ~  . + ",bspl, "*",interactionTerm, ")", sep="")
    currentModel<-eval(parse(text=test)) 
  }
  
  
  tempFit <- get.measure_2d(fitnessMeasure, currentFit, currentModel,data, dists,aR,radii,radiusIndices)$fitStat
  if(tempFit <= (currentFit+10)){
    models[[length(models)+1]] = list(aR,radiusIndices, radii, tempFit)
  }
  modelinprogress<<-currentModel
  #print("ooooooooooooooooooooooooooooooooooooooo")
  #print("Model fitted...")
  #print("ooooooooooooooooooooooooooooooooooooooo")
  return(list(currentModel=currentModel,models=models))
  
}