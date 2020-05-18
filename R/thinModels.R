#' function to thin the number of models
#' 
#' @author Cameron Walker, Department of Enginering Science, University of Auckland.
#' 
#' @export
#' 

thinModels<- function(models){
  
  if(length(models)>0){
    bic = vector(length=length(models))
    
    for(i in 1:length(models)){
      bic[i] = models[[i]][[4]]
    }
    
    id = which(bic<=(min(bic)+10))
    tempModels = models[id]  
  }else{
    tempModels = models
  }
  return(models = tempModels)
}


thinModels.hr<- function(models, submodels){
  
  if(length(models)>0){
    bic = vector(length=length(models))
    
    for(i in 1:length(models)){
      bic[i] = models[[i]][[4]] + submodels[[i]][[1]]
    }
    
    id = which(bic<=(min(bic)+10))
    tempModels = models[id] 
    tempSub = submodels[id]
  }else{
    tempModels = models
    tempSub = submodels
  }
  return(list(models = tempModels, submodels=tempSub))
}